#include <cpl.h>

cpl_parameterlist* fillParList(void);
void loadData(cpl_table** spec, cpl_size nspec, char* datapath);
cpl_propertylist* loadPropList(char* datapath);
cpl_table* fillMolecules(void);
cpl_table* fillExcludedPixels(void);
