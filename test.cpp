#include "solver.hpp"
#include "odeset.hpp"

int main()
{
  auto solver = solver_t<odeset_t>(odeset_t());
  while (solver.t < 200)
  {
    solver.step();
  }
}
