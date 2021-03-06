<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark pour évaluer le calcul des Cqs et la maj du vecteur à partir du tenseur et la maj multi-env du tenseur</title>
    <timeloop>CqsVectorTensorLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>Nbenv</variable>
      <variable>VolumeVisu</variable>
      <variable>Volume</variable>
      <variable>FracVol</variable>
      <variable>FracVolVisu</variable>
<!--       <variable>Tensor</variable> -->
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
        <nsd>2 2 1</nsd>
        <origine>0. 0. 0.</origine>
        <lx nx="30" prx="1.0">1.</lx>
        <ly ny="30" pry="1.0">1.</ly>
        <lz nz="30" pry="1.0">1.</lz>
      </cartesian>
    </meshgenerator>
  </mesh>

  <!-- Configuration du module GeomEnv -->
  <geom-env>
    <visu-volume>true</visu-volume>
    <visu-frac-vol>true</visu-frac-vol>
    <geom-scene>env5m3</geom-scene>
  </geom-env>
</case>
