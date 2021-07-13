#include <arcane/launcher/ArcaneLauncher.h>

using namespace Arcane;

int
main(int argc,char* argv[])
{
//#define DEPRECATED
#ifdef DEPRECATED
  auto& app_info = ArcaneLauncher::applicationInfo();
  app_info.setCommandLineArguments(CommandLineArguments(&argc,&argv));
  app_info.setCodeName("Pattern4GPU");
  return ArcaneLauncher::run();
#else
  ArcaneLauncher::init(CommandLineArguments(&argc,&argv));
  ApplicationBuildInfo& app_build_info = ArcaneLauncher::applicationBuildInfo();
  app_build_info.setCodeName("Pattern4GPU");
  return ArcaneLauncher::run();
#endif
}
