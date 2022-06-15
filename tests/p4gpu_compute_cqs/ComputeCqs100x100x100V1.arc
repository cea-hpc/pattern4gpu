<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark pour évaluer le calcul des Cqs sur allCells() et la maj du vecteur sur active_cells</title>
    <timeloop>ComputeCqsAndVectorLoop</timeloop>
  </arcane>

<!--   <arcane-post-processing> -->
<!--     <output-period>1</output-period> -->
<!--     <output> -->
<!--       <variable>Nbenv</variable> -->
<!--       <variable>VolumeVisu</variable> -->
<!--       <variable>Volume</variable> -->
<!--     </output> -->
<!--     <format> -->
<!--       <binary-file>false</binary-file> -->
<!--     </format> -->
<!--   </arcane-post-processing> -->

  <!-- ***************************************************************** -->
  <!--Definition du maillage cartesien -->
  <mesh nb-ghostlayer="3" ghostlayer-builder-version="3">
    <meshgenerator>
      <cartesian>
        <nsd>2 2 1</nsd>
        <origine>0. 0. 0.</origine>
        <lx nx="100" prx="1.0">1.</lx>
        <ly ny="100" pry="1.0">1.</ly>
        <lz nz="100" pry="1.0">1.</lz>
      </cartesian>
    </meshgenerator>
  </mesh>

  <!-- Configuration du module GeomEnv -->
  <geom-env>
    <visu-volume>false</visu-volume>
    <geom-scene>env5m3</geom-scene>
  </geom-env>

  <!-- Configuration du service AccEnvDefault -->
  <acc-env-default>
    <acc-mem-advise>true</acc-mem-advise>
    <device-affinity>node_rank</device-affinity>
    <!-- <heterog-partition>none</heterog-partition> -->
  </acc-env-default>

  <!-- Configuration du module Pattern4GPU -->
  <pattern4-g-p-u>

    <init-cqs-version>arcgpu_v1</init-cqs-version>
    <init-node-vector-version>arcgpu_v1</init-node-vector-version>
    <init-node-coord-bis-version>arcgpu_v1</init-node-coord-bis-version>
    <init-cell-arr12-version>arcgpu_v1</init-cell-arr12-version>
    <!-- <compute-cqs-vector-version>ori</compute-cqs-vector-version> -->
    <compute-cqs-vector-version>arcgpu_v1</compute-cqs-vector-version>
  </pattern4-g-p-u>
</case>
