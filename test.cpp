#include <algorithm>
#include <limits>

#include "src/solver.hpp"
#include "src/odeset.hpp"

int main()
{
  auto odeset = odeset_t();
  auto solver = solver_t<odeset_t>(odeset);

  auto t_last = solver.t, dtmin = std::numeric_limits<double>::infinity();
  while (solver.t < 200)
  {
    solver.step();
    auto RH = odeset.RH(solver.state, solver.t);
std::cout 
  << solver.t
  << " " 
  << RH
  << std::endl;
    if (RH > 1) dtmin = std::min(dtmin, solver.t - t_last);
    t_last = solver.t;
  }
std::cerr << dtmin <<  std::endl;
}
