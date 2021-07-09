<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Calculer par direction les volumes à travers les faces</title>
    <timeloop>CartVolLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>DirVol2Left</variable>
      <variable>DirVol2Right</variable>
      <variable>FaceVelocityRight</variable>
      <variable>DirDefCoordLeft</variable>
      <variable>DirTransAreaLeft</variable>
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
  </arcane-post-processing>

  <!-- ***************************************************************** -->
  <!--Definition du maillage cartesien -->
  <mesh nb-ghostlayer="3" ghostlayer-builder-version="3">
    <meshgenerator>
      <cartesian>
        <nsd>4 3 2</nsd>
        <origine>0. 0. 0.</origine>
        <lx nx="10" prx="1.0">1.</lx>
        <ly ny="8" pry="1.0">1.</ly>
        <lz nz="5" pry="1.0">1.</lz>
      </cartesian>
    </meshgenerator>
  </mesh>

</case>
