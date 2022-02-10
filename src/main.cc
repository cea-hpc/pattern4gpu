#include <arcane/launcher/ArcaneLauncher.h>
#include <arcane/utils/Exception.h>

using namespace Arcane;

// Défini dans accenv/MsgPassInit.cc
void msg_pass_init(int*, char***);

int
main(int argc,char* argv[])
{
  try {
    msg_pass_init(&argc, &argv);

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
}
