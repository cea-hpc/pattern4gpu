#include <arcane/launcher/ArcaneLauncher.h>
#include <arcane/utils/Exception.h>

using namespace Arcane;

int
main(int argc,char* argv[])
{
  try {
    ArcaneLauncher::init(CommandLineArguments(&argc,&argv));
    ApplicationBuildInfo& app_build_info = ArcaneLauncher::applicationBuildInfo();
    //app_build_info.setMessagePassingService("SequentialParallelSuperMng");
    app_build_info.setCodeName("Pattern4GPU");
    return ArcaneLauncher::run();
  }
  catch(const Arcane::Exception& ex){
    std::cerr << "EXCEPTION: " << ex << "\n";
    return 1;
  }
  return 0;  
}
