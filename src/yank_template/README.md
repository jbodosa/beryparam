## Steps to generate the alchem setup

<ol>
  <li>Generate the moleucle/system with counter-ions using packmol > add_ions.py</li>
  <li>Generate the crd file for this ion box</li>
  <li>Use the crd file to generate the water and vacuum boxes using CHARMM</li>

```
$PACKMOL <ion_box.inp >ion_box.out

$CHARMM <ion_wat.inp >ion_wat.out
$CHARMM <ion_vac.inp >ion_vac.out
```

## Steps
1. Recenter : recenter and convert the molecule/system pdb to CHARMM pdb, psf, crd ([charmm const](https://www.charmm-gui.org/charmmdoc/subst.html))
2. calculate charge of molecule/system and add charge accordingly at a specified distance, output pdb, psf, crd
3. Add restrain to the counter-ions and create xml file of restrained system
4. Probably also add a cleanme script and a testme script


## NOTES
- Problems adding water using charmm ff  [link](https://github.com/openmm/openmm/issues/3566)
- MDAnalysis can ead and write crd but it might require segid in pdb or the psf. And it does not write the charmm-extended format.

