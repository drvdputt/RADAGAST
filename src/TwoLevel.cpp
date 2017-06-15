#include "TwoLevel.h"
#include "Constants.h"
#include "Sanity.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "TwoLevelHardcoded.h"
#include "global.h"
#include <algorithm>
#include <exception>
#include <vector>

TwoLevel::TwoLevel()
                : NLevel(new TwoLevelHardcoded())
{
}
