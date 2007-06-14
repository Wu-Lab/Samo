
#include "Samo.h"

#include "MemLeak.h"


int main(int argc, char *argv[])
{
	EnableMemLeakCheck();
//	_CrtSetBreakAlloc(1682);

	try {
		Samo samo(argc, argv);
		samo.run();
	}
    catch (exception &e)
    {
		Logger::error("%s", e.what());
        return 1;
	}
	return 0;
}
