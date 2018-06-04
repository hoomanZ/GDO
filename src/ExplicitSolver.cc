#include <Solver.h>

using namespace FiniteElement;

Solver::Solver()
    :d_time(1.0), d_del_t(0.001), d_box(new Box())
{
  setSolverType(Solver::solverType::Explicit);
}

Solver::Solver(solverType type, double time, double delT, BoxSP box)
    :d_time(time), d_del_t(delT), d_box(box)
{
//setBox(box);
  setSolverType(type);
  if (type == Explicit)
     explicitSolver();
  else 
     implicitSolver();
}

Solver::~Solver()
{
}


void
Solver::explicitSolver()
{
  d_box->makeGlobalMatrices();
  std::vector<double> extForce = d_box->getExternalForce(); 
  std::vector<double> oldDisp = d_box->getOldGlobDisp();
  std::vector<double> disp = d_box->getGlobDisp();
  std::vector<double> newDisp; // = d_box->getNewGlobDisp();

  std::vector<std::vector<double> > stiffMatrix = d_box->getGlobalStiffness();
  std::vector<std::vector<double> > massLumpedMatrix = d_box->getGlobalLumpedMass();
  int size = d_box->getGlobalStiffness().size();
  int iter_num = 1;
  while (iter_num*d_del_t < d_time)
  {
    newDisp.clear();
    newDisplacement(stiffMatrix, massLumpedMatrix, extForce, oldDisp, disp, newDisp);
    std::cout << "iter_num= " << iter_num;
    for (int ll = 0; ll < newDisp.size(); ll++)
    {
      if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
      std::cout << "(" << ll << ")= " << newDisp[ll] << "  ";
    }
    std::cout << std::endl;
//    d_box->getOldGlobDisp().clear();
//    d_box->initialVector(d_box->getOldGlobDisp(), size);
//    d_box->getGlobDisp().clear();
//    d_box->initialVector(d_box->getGlobDisp(), size);
    d_box->setOldGlobDisp(disp);
    d_box->setGlobDisp(newDisp);
    oldDisp = d_box->getOldGlobDisp();
    disp = d_box->getGlobDisp();
    iter_num++;
  }
}

void
Solver::newDisplacement(const std::vector<std::vector<double> >& stiffMatrix,
                                const std::vector<std::vector<double> >& massLumpedMatrix,
                                const std::vector<double>& extForce, const std::vector<double>& oldDisp,
                                const std::vector<double>& disp, std::vector<double>& newDisp)
{
  int size = d_box->getGlobalStiffness().size();
  double sum = 0.0;
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      sum += stiffMatrix[i][j]*disp[j];
    }
    newDisp.push_back((d_del_t*d_del_t)*(1/massLumpedMatrix[i][i])*(extForce[i]-sum)+2*disp[i]-oldDisp[i]);
  }
}

void 
Solver::implicitSolver()
{
}









