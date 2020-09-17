#include <ibtk/AppInitializer.h>
#include <ibtk/FEDataManager.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/app_namespaces.h>

#include <tbox/Logger.h>

#include "../tests.h"

int
main(int argc, char** argv)
{
    using namespace IBTK;
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Since this is a test we do not want to print file names or line numbers
    // to output files:
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    SubdomainToPatchLevelTranslation trans(10, input_db->getDatabase("trans"));

    for (unsigned int i = 0; i < 20; ++i)
    {
        plog << "subdomain id: " << i << " level: " << trans[i] << '\n';
    }

    for (unsigned int i = 1000; i < 1010; ++i)
    {
        plog << "subdomain id: " << i << " level: " << trans[i] << '\n';
    }
}
