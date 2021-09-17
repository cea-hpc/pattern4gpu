<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark pour évaluer le calcul des Cqs sur allCells() et la maj du vecteur sur active_cells</title>
    <timeloop>ComputeCqsAndVectorLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>Nbenv</variable>
      <variable>VolumeVisu</variable>
      <variable>Volume</variable>
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
        <lx nx="10" prx="1.0">1.</lx>
        <ly ny="10" pry="1.0">1.</ly>
        <lz nz="10" pry="1.0">1.</lz>
      </cartesian>
    </meshgenerator>
  </mesh>

  <!-- Configuration du module GeomEnv -->
  <geom-env>
    <visu-volume>true</visu-volume>
    <geom-scene>env5m3</geom-scene>
  </geom-env>

  <!-- Configuration du module Pattern4GPU -->
  <pattern4-g-p-u>
    <acc-mem-advise>true</acc-mem-advise>

    <init-cell-arr12-version>arcgpu_v1</init-cell-arr12-version>
    <!-- <compute-cqs-vector-version>ori</compute-cqs-vector-version> -->
    <compute-cqs-vector-version>arcgpu_v1</compute-cqs-vector-version>
  </pattern4-g-p-u>
</case>
