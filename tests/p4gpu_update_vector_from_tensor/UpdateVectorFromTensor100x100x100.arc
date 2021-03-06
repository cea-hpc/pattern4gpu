<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark pour évaluer maj grandeurs tenseurs multi-env</title>
    <timeloop>UpdateVectorFromTensorLoop</timeloop>
  </arcane>

<!--   <arcane-post-processing> -->
<!--     <output-period>1</output-period> -->
<!--     <output> -->
<!--       <variable>Nbenv</variable> -->
<!--       <variable>VolumeVisu</variable> -->
<!--       <variable>Tensor</variable> -->
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

<!--   <arcane-checkpoint> -->
<!--     <period>0</period> -->
    <!-- Mettre '0' si on souhaite ne pas faire de protections a la fin du calcul -->
<!--     <do-dump-at-end>0</do-dump-at-end> -->
<!--     <checkpoint-service name="ArcaneBasic2CheckpointWriter" /> -->
<!--   </arcane-checkpoint> -->

  <!-- Configuration du module GeomEnv -->
  <geom-env>
    <visu-volume>false</visu-volume>
    <geom-scene>env5m3</geom-scene>
  </geom-env>
  
  <!-- Configuration du module Pattern4GPU -->
  <pattern4-g-p-u>
    <init-node-vector-version>ori</init-node-vector-version>
    <!-- <init-node-vector-version>mt</init-node-vector-version> -->
    <init-cqs1-version>arcgpu_v1</init-cqs1-version>
    <update-vector-from-tensor-version>ori</update-vector-from-tensor-version>
    <!-- <update-vector-from-tensor-version>mt</update-vector-from-tensor-version> -->
  </pattern4-g-p-u>
</case>
