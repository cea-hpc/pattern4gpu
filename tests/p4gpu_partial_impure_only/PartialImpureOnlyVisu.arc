<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark la maj de valeurs partielles sur les mailles mixtes (impures)</title>
    <timeloop>PartialImpureOnlyLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>Nbenv</variable>
      <variable>FracVolVisu</variable>
      <variable>FracVol</variable>
      <variable>MEnvVar1Visu</variable>
      <variable>MEnvVar1</variable>
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

  <!-- Configuration du module Pattern4GPU -->
  <pattern4-g-p-u>
    <acc-mem-advise>true</acc-mem-advise>

    <visu-m-env-var>true</visu-m-env-var>
    <partial-impure-only-version>ori</partial-impure-only-version>
  </pattern4-g-p-u>
</case>
