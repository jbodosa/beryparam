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
test
