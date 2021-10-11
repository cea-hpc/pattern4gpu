<?xml version='1.0'?>
<case codeversion="1.0" codename="Pattern4GPU" xml:lang="en">
  <arcane>
    <title>Benchmark la maj de valeurs partielles sur les mailles mixtes (impures)</title>
    <timeloop>PartialImpureOnlyLoop</timeloop>
  </arcane>

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
    <visu-frac-vol>false</visu-frac-vol>
    <!-- <geom-scene>env5m3</geom-scene> -->
    <geom-scene>nestNdiams</geom-scene>
    <nested-ndiams>9</nested-ndiams>
  </geom-env>

  <!-- Configuration du module Pattern4GPU -->
  <pattern4-g-p-u>
    <acc-mem-advise>true</acc-mem-advise>

    <visu-m-env-var>false</visu-m-env-var>
    <init-menv-var-version>ori</init-menv-var-version>
    <!-- <init-menv-var-version>arcgpu_v1</init-menv-var-version> -->
    <partial-impure-only-version>ori</partial-impure-only-version>
    <!-- <partial-impure-only-version>arcgpu_v1</partial-impure-only-version> -->
  </pattern4-g-p-u>
</case>
