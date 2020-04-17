#include "CollisionParameters.hpp"
#include "DoctestUtils.hpp"
#include "H2Data.hpp"
#include "SimpleColumnFile.hpp"
#include "SpeciesIndex.hpp"
#include "Testing.hpp"

using namespace GasModule;

TEST_CASE("Collision coefficients")
{
    // this will load and write data to try and recreate figure 2 from Lee et al. 2008
    H2Data h2data;
    SpeciesIndex spindex = SpeciesIndex::makeDefault();
    SpeciesVector sv(&spindex);

    // c = rate coefficient (q) * collision partner density (nH2). By setting nH2 to 1, c should
    // equal q I think.
    sv.setNH2(1);

    // range of temperatures for the plot
    Array tv = Testing::generateGeometricGridv(200, 2, 10000);

    std::vector<std::string> columns = {"temperature", "J2-J0", "J4-J2", "J6-J4", "J8-J6"};
    std::vector<std::array<int, 2>> jtransitions = {{2, 0}, {4, 2}, {6, 4}, {8, 6}};

    // pure para
    if (DoctestUtils::allowFileOutput)
    {
        OutColumnFile out("lee08rate-para-para", columns, 6);
        for (double t : tv)
        {
            CollisionParameters cp(t, sv, 0);
            EMatrix cvv = h2data.cvv(cp);

            // gather data for line
            std::vector<double> row;
            row.reserve(columns.size() + 1);
            row.push_back(t);

            for (auto jul : jtransitions)
            {
                int upper = h2data.indexFind(H2Data::ElectronicState::X, jul[0], 0);
                int lower = h2data.indexFind(H2Data::ElectronicState::X, jul[1], 0);
                row.push_back(cvv(upper, lower));
            }
            out.writeLine(row);
        }
    }
}
