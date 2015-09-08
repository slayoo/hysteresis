#include "solver.hpp"
#include "odeset.hpp"

int main()
{
  auto solver = solver_t<odeset_t>(odeset_t());
  solver.step();
}
