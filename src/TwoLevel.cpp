#include "TwoLevel.h"
#include "Constants.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "TwoLevelHardcoded.h"
#include "global.h"
#include <algorithm>
#include <exception>
#include <vector>
#include "Error.h"

TwoLevel::TwoLevel()
                : NLevel(new TwoLevelHardcoded())
{
}
