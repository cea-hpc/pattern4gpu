<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark </title>
    <timeloop>PartialAndMean4Loop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>Nbenv</variable>
      <variable>FracVolVisu</variable>
      <variable>FracVol</variable>
      <variable>MEnvVar1Visu</variable>
      <variable>MEnvVar1</variable>
      <variable>MEnvVar2Visu</variable>
      <variable>MEnvVar2</variable>
      <variable>MEnvVar3Visu</variable>
      <variable>MEnvVar3</variable>
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
    <visu-frac-vol>true</visu-frac-vol>
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

    <visu-m-env-var>true</visu-m-env-var>
    <init-menv-var-version>ori</init-menv-var-version>
    <partial-and-mean4-version>ori</partial-and-mean4-version>
  </pattern4-g-p-u>
</case>
