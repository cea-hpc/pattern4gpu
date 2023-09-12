<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Calculer par direction les volumes à travers les faces</title>
    <timeloop>CartVolLoop</timeloop>
  </arcane>

  <!-- ***************************************************************** -->
  <!--Definition du maillage cartesien -->
  <mesh nb-ghostlayer="3" ghostlayer-builder-version="3">
    <meshgenerator>
      <cartesian>
        <nsd>2 2 2</nsd>
        <origine>0. 0. 0.</origine>
        <lx nx="100" prx="1.0">1.</lx>
        <ly ny="100" pry="1.0">1.</ly>
        <lz nz="100" pry="1.0">1.</lz>
      </cartesian>
    </meshgenerator>
  </mesh>

  <acc-env-default>
    <acc-mem-advise>true</acc-mem-advise>
  </acc-env-default>

  <!-- Configuration du module Pattern4GPU -->
  <pattern4-g-p-u>

    <visu-m-env-var>false</visu-m-env-var>
    <!-- <compute-vol-version>ori</compute-vol-version> -->
    <!-- <compute-vol-version>arcgpu_v1</compute-vol-version> -->
    <compute-vol-version>arcgpu_v2</compute-vol-version>
  </pattern4-g-p-u>

</case>
