#include <arcane/launcher/ArcaneLauncher.h>
#include <arcane/utils/Exception.h>
#include <mpi.h>

using namespace Arcane;

int
main(int argc,char* argv[])
{
  try {
//#define USE_THREAD_MULTIPLE
#ifdef USE_THREAD_MULTIPLE
    // On impose le niveau MPI_THREAD_MULTIPLE
    Integer provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ARCANE_ASSERT(provided==MPI_THREAD_MULTIPLE, 
        ("Impossible d'utiliser le niveau MPI_THREAD_MULTIPLE"));
    Integer rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0) {
      if (provided==MPI_THREAD_SERIALIZED)
        std::cout << "MPI_THREAD_SERIALIZED" << std::endl;
      else if (provided==MPI_THREAD_MULTIPLE)
        std::cout << "MPI_THREAD_MULTIPLE" << std::endl;
      else
        std::cout << "MPI_THREAD_{SINGLE|FUNNELED}" << std::endl;
    }
    // Fin init MPI en dur
#endif

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
