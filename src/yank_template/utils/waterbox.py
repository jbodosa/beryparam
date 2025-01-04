import os
from datetime import datetime, timezone
import logging

import re

logger = logging.getLogger(__name__)

ION_LIST = ["SOD", "POT", "CAL","CLA"]

### 1. Write the CHARMM inp to generate the waterbox
### 2. Use the waterbox to solvate the non-charged system and remove overlapping water
### 3. Use the waterbox to add ions and solvate solvate the charged system and remove overlapping water

def write_waterbox_inp(xdim, ydim, zdim):
    """
        Write the CHARMM inp to generate waterbox of ceratin size
    """
    output_file = "generate_waterbox.inp"

    with open(output_file, 'w') as inp:
        # Get username
        user = os.getenv('USERNAME') or os.getenv('USER')
        # Get current UTC date and time
        current_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")

        inp.write(f"* GENERATED BY SYSGEN\n")
        inp.write(f"* THIS IS TO GENRATE A WATERBOX OF GIVEN SIZE\n")
        inp.write(f"* Made by {user} on {current_utc}\n")
        inp.write(f"* \n")
        inp.write(f"\n")

        inp.write(f"DIMENS CHSIZE 5000000 MAXRES 3000000\n")

        inp.write(f"! Read topology and parameter files\n")
        inp.write(f"stream toppar.str\n")

        inp.write(f"calc Xinit = {xdim} !from user' input\n")
        inp.write(f"calc Yinit = {ydim}\n")
        inp.write(f"calc Zinit = {zdim}\n")
        inp.write(f"calc Rinit = 0\n")

        inp.write(f"calc xcen = 0.0\n")
        inp.write(f"calc ycen = 0.0\n")
        inp.write(f"calc zcen = 0.0\n")
        inp.write(f"! parameters for water box & image (crystal)\n")

        inp.write(f"calc BoxSizeX = @Xinit\n")
        inp.write(f"calc BoxSizeY = @Yinit\n")
        inp.write(f"calc BoxSizeZ = @Zinit\n")

        inp.write(f"set XTLtype = ORTHorhombic\n")
        inp.write(f"if BoxSizeX .eq. @BoxSizeY if BoxSizeX .ne. @BoxSizeZ set XTLtype = TETRagonal\n")
        inp.write(f"if BoxSizeX .eq. @BoxSizeY if BoxSizeX .eq. @BoxSizeZ set XTLtype = CUBIc\n")
        inp.write(f"set A = @BoxSizeX\n")
        inp.write(f"set B = @BoxSizeY\n")
        inp.write(f"set C = @BoxSizeZ\n")
        inp.write(f"set Alpha = 90.0\n")
        inp.write(f"set Beta  = 90.0\n")
        inp.write(f"set Gamma = 90.0\n")

        inp.write(f"! a pre-equilibrated water cubic box with L=18.8560\n")
        inp.write(f"set L 18.8560\n")

        inp.write(f"! number of boxes along XYZ-directions\n")
        inp.write(f"calc Xnum = int(@BoxSizeX/@L) + 1\n")
        inp.write(f"calc Ynum = int(@BoxSizeY/@L) + 1\n")
        inp.write(f"calc Znum = int(@BoxSizeZ/@L) + 1\n")

        inp.write(f"! base unit of water box\n")
        inp.write(f"read sequence TIP3 216\n")
        inp.write(f"generate W000 setup noangle nodihedral\n")

        inp.write(f"open read unit 10 card name tip216.crd\n")
        inp.write(f"read coor unit 10 card\n")
        inp.write(f"close unit 10\n")


        inp.write(f"coor stat sele type OH2 end\n")
        inp.write(f"calc Lhalf = @L / 2.0\n")
        inp.write(f"coor trans xdir 1.0 dist @Lhalf\n")
        inp.write(f"coor trans ydir 1.0 dist @Lhalf\n")
        inp.write(f"coor trans zdir 1.0 dist @Lhalf\n")
        inp.write(f"coor stat sele type OH2 end\n")

        inp.write(f"! planar water box unit (XY)\n")

        inp.write(f"set J2  1\n")
        inp.write(f"label DO_2\n")
        inp.write(f"    set J1  1\n")
        inp.write(f"    label DO_1\n")

        inp.write(f"    calc wsegid = ( @J2 - 1 ) * @Xnum + @J1\n")

        inp.write(f"    read sequence TIP3 216\n")
        inp.write(f"    generate W@wsegid setup noangle nodihedral\n")

        inp.write(f"    coor duplicate select segid W000 end select segid W@wsegid end\n")

        inp.write(f"    calc X = @L * ( @J1 - 1 )  ! mult X by @J1\n")
        inp.write(f"    calc Y = @L * ( @J2 - 1 )  ! mult Y by @J2\n")

        inp.write(f"    coor trans xdir @X ydir @Y select segid W@wsegid end\n")

        inp.write(f"    incr J1 by 1\n")
        inp.write(f"    if J1 le @Xnum goto DO_1\n")
        inp.write(f"incr J2 by 1\n")
        inp.write(f"if J2 le @Ynum goto DO_2\n")

        inp.write(f"define junk sele .byres. ( ( type OH2 ) .and. -\n")
        inp.write(f"                           ( prop X .gt. @BoxSizeX .or. -\n")
        inp.write(f"                             prop Y .gt. @BoxSizeY ) ) end\n")
        inp.write(f"if ?nsel .gt. 0 delete atoms sele junk end\n")

        inp.write(f"define solv sele .byres. type OH2 end\n")
        inp.write(f"if ?nsel .eq. 0 stop ! ABNORMAL TERMINATION: Too small box size\n")

        inp.write(f"delete atom sele segid W000 end\n")

        inp.write(f"open write unit 10 card name water_tmp.crd\n")
        inp.write(f"write coor unit 10 card\n")

        inp.write(f"delete atom sele all end\n")

        inp.write(f"! generate water box by stacking planar water boxes along Z\n")

        inp.write(f"set J3  1\n")
        inp.write(f"label DO_3\n")

        inp.write(f"    open read card unit 10 name water_tmp.crd\n")
        inp.write(f"    read sequence coor card unit 10\n")
        inp.write(f"    generate Wz@J3 setup warn noangle nodihedral\n")

        inp.write(f"    open read unit 10 card name water_tmp.crd\n")
        inp.write(f"    read coor unit 10 card append\n")

        inp.write(f"    calc Z = @L * ( @J3 - 1 )  ! mult Y by @J3\n")
        inp.write(f"    coor trans zdir @Z select segid Wz@J3 end\n")

        inp.write(f"incr J3 by 1\n")
        inp.write(f"if J3 le @Znum goto DO_3\n")

        inp.write(f"define junk sele .byres. ( ( type OH2 ) .and. -\n")
        inp.write(f"                           ( prop Z .gt. @BoxSizeZ ) ) end\n")
        inp.write(f"if ?nsel .gt. 0 delete atoms sele junk end\n")

        inp.write(f"coor stat sele type OH2 end\n")
        inp.write(f"coor orient norotation\n")
        inp.write(f"coor stat sele type OH2 end\n")

        inp.write(f"!\n")
        inp.write(f"!Shaping the box\n")
        inp.write(f"!\n")

        inp.write(f"COOR CONVERT ALIGNED SYMMETRIC @A @B @C @alpha @beta @gamma\n")
        inp.write(f"coor copy comp\n")

        inp.write(f"CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma\n")
        inp.write(f"CRYSTAL BUILD NOPER 0 CUTOFF 2.0\n")

        inp.write(f"!Image centering by residue\n")
        inp.write(f"IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen\n")

        inp.write(f"nbond ctonnb 2.0 ctofnb 3.0 cutnb 3.0 cutim 3.0 wmin 0.001\n")
        inp.write(f"CRYSTAL FREE\n")

        inp.write(f"coor diff comp\n")
        inp.write(f"define junk sele .byres. ( ( prop Xcomp .ne. 0 ) .or. -\n")
        inp.write(f"                           ( prop Ycomp .ne. 0 ) .or. -\n")
        inp.write(f"                           ( prop Zcomp .ne. 0 ) ) end\n")
        inp.write(f"if ?nsel .gt. 0 delete atoms sele junk end\n")

        inp.write(f"COOR CONVERT SYMMETRIC ALIGNED @A @B @C @alpha @beta @gamma\n")


        inp.write(f"coor stat sele type OH2 end\n")
        inp.write(f"set nwater ?nsel\n")

        inp.write(f"open write unit 10 card name waterbox.crd\n")
        inp.write(f"write coor unit 10 card\n")
        inp.write(f"* Equilibrated water\n")
        inp.write(f"*\n")

        inp.write(f"open write card unit 2 name waterbox.pdb\n")
        inp.write(f"write coor pdb  unit 2\n")
        inp.write(f"* Equilibrated water\n")
        inp.write(f"*\n")

        inp.write(f"open write unit 90 card name waterbox.str\n")
        inp.write(f"write title unit 90\n")
        inp.write(f"* read sequence TIP3 @nwater\n")
        inp.write(f"* generate SOLV setup noangle nodihedral\n")
        inp.write(f"*\n")

        inp.write(f"open  write unit 90 card name waterbox.prm\n")
        inp.write(f"write title unit 90\n")
        inp.write(f"* set BoxType  = rect\n")
        inp.write(f"* set XTLtype  = @XTLtype\n")
        inp.write(f"* set A = @A\n")
        inp.write(f"* set B = @B\n")
        inp.write(f"* set C = @C\n")
        inp.write(f"* set Alpha = @Alpha\n")
        inp.write(f"* set Beta  = @Beta\n")
        inp.write(f"* set Gamma = @Gamma\n")
        inp.write(f"* set xcen = @xcen\n")
        inp.write(f"* set ycen = @ycen\n")
        inp.write(f"* set zcen = @zcen\n")
        inp.write(f"*\n")

        inp.write(f"stop\n")
    return()

def solvate_system(mol_segid_list, mol_n_list):
    """
    Solvate the un-charged system.
    """
    output_file = "solvate_system.inp"

    with open(output_file, 'w') as inp:
        # Get username
        user = os.getenv('USERNAME') or os.getenv('USER')
        # Get current UTC date and time
        current_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")

        inp.write(f"* GENERATED BY SYSGEN\n")
        inp.write(f"* THIS IS TO GENRATE A WATERBOX OF GIVEN SIZE\n")
        inp.write(f"* Made by {user} on {current_utc}\n")
        inp.write(f"* \n")
        inp.write(f"\n")

        inp.write(f"DIMENS CHSIZE 5000000 MAXRES 3000000\n")

        inp.write(f"! Read topology and parameter files\n")
        inp.write(f"stream toppar.str\n")

        mol_segid_sel = "( "
        for mol_segid, mol_n in  zip(mol_segid_list, mol_n_list):
            mol_segid_sel+=f" segid {mol_segid} .or."
            if mol_segid not in ION_LIST:
                # If not ion
                inp.write(f"! read the sequence\n")
                inp.write(f"read sequence {mol_segid} {mol_n}\n")
                inp.write(f"generate {mol_segid}\n")
                inp.write(f"\n")
            else:
                # If it is ion
                inp.write(f"! read the sequence\n")
                inp.write(f"read sequence {mol_segid} {mol_n}\n")
                inp.write(f"generate {mol_segid} noangle nodihedral\n")
                inp.write(f"\n")
        #        mol_segid_sel.removesuffix(".or.")
        mol_segid_sel = re.sub('\.or.$', '', mol_segid_sel)
        mol_segid_sel+=" )"

        inp.write(f"! Read the coordinates of the non-charged molecule\n")
        inp.write(f"open read card unit 31 name recenter.crd\n")
        inp.write(f"read coor card unit 31 \n")
        inp.write(f"close unit 31\n")

        inp.write(f"! Read water box parameters and coordinates\n")
        inp.write(f"stream waterbox.str\n")

        inp.write(f"open read card unit 30 name waterbox.crd\n")
        inp.write(f"read coor card unit 30 append\n")
        inp.write(f"close unit 30\n")

        inp.write(f"! Remove water molecules overlapping with the non-charged molecule\n")
        inp.write(f"define NONCHG sele .not. hydrogen .and. {mol_segid_sel} end\n")
        inp.write(f"if ?nsel .gt. 0 then\n")
        inp.write(f"    define target sele .byres. ( ( type OH2 .and. segid SOLV ) .and. ( NONCHG .around. 2.8 ) ) end\n")
        inp.write(f"    if ?nsel .gt. 0 delete atom sele target end\n")
        inp.write(f"endif\n")

        inp.write(f"! Check and count the remaining water molecules\n")
        inp.write(f"coor stat sele type OH2 .and. segid SOLV end\n")
        inp.write(f"set nwater ?nsel\n")

        inp.write(f"! Write coordinates and system information\n")
        inp.write(f"open write unit 10 card name solvated_sys.psf\n")
        inp.write(f"write psf unit 10 card\n")

        inp.write(f"open write card unit 10 name solvated_sys.pdb\n")
        inp.write(f"write coor pdb unit 10\n")

        inp.write(f"open write unit 10 card name solvated_sys.crd\n")
        inp.write(f"write coor unit 10 card\n")

        inp.write(f"open write unit 90 card name solvated_sys.str\n")
        inp.write(f"write title unit 90\n")
        inp.write(f"** Solvated system settings\n")
        inp.write(f"**\n")
        inp.write(f"* read sequence TIP3 @nwater\n")
        inp.write(f"* generate SOLV setup noangle nodihedral\n")
        inp.write(f"*\n")

        inp.write(f"stop\n")

    return()
