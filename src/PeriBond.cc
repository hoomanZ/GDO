#include <PeriBond.h>
#include <PeriMaterialPoint.h>

using namespace FiniteElement;


PeriBond::PeriBond()
  : d_center_point(new PeriMaterialPoint()), 
    d_second_point(new PeriMaterialPoint()),
    d_influence_function(1.0),
    d_broken_bond(false)
{
  Vector3D vec(1.0, 1.0, 1.0);
  setPairwiseForce(vec);
}


PeriBond::PeriBond(PeriMaterialPointP center, PeriMaterialPointP second, bool broken)
  : d_center_point(center), d_second_point(second), 
    d_influence_function(1.0), d_broken_bond(false)
{
  Vector3D vec(1.0, 1.0, 1.0);
  setPairwiseForce(vec);
}


void
PeriBond::clearBond()
{
  PeriMaterialPointP center(new PeriMaterialPoint());
  PeriMaterialPointP second(new PeriMaterialPoint());
  setCenterPoint(center);  
  setSecondPoint(second);
  setBrokenBond(false);
  Vector3D vec(1.0, 1.0, 1.0);
  setPairwiseForce(vec);
  setInfluenceFunction(1.0);

}



double
PeriBond::lengthOfSai()
{
  return (d_second_point->getPosOld() - d_center_point->getPosOld()).length();
}



void
PeriBond::calculateVolume(double horizon, double rj, double delX)
{
  Vector3D sai = getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld();
  double distance = sai.length();  
  double volume = delX*delX*delX;
  double coefficient = 0.5 + (horizon - distance)/(2*rj);
  if (((horizon - rj) <= distance) && (distance <= horizon))
  {
    volume *= coefficient;
  }
  else if (distance <= (horizon - rj))
  {
    volume *= 1;
  }
  else 
  {
    volume = 0.0;
  }

  setVolume(volume);

}




