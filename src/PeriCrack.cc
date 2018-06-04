#include <PeriCrack.h>

#include<PeriMaterialPointPArray.h>
#include <PeriMaterialPointP.h>
#include <PeriMaterialPoint.h>

#include<PeriBondPArray.h>
#include<PeriBondP.h>
#include<PeriBond.h>

#include <StateBondPArray.h>
#include <StateBondP.h>
#include <StateBond.h>

using namespace FiniteElement;


PeriCrack::PeriCrack()
  :d_min_point(0.0, 0.0, 0.0), d_max_point(0.1, 0.1, 0.1)
{


}


PeriCrack::PeriCrack(Vector3D crackMinPoint, Vector3D crackMaxPoint)
{
  setMinPoint(crackMinPoint);
  setMaxPoint(crackMaxPoint);
}



void
PeriCrack::createInitialCrack(PeriMaterialPointPArray& pointsArray)
{
  Vector3D zero(0.0, 0.0, 0.0);
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_point = pointsArray[i];
    PeriBondPArray bondsArray = cur_point->getHorizonBonds();
    for (int j = 0; j < bondsArray.size(); j++)
    {
      PeriBondP cur_bond = bondsArray[j];
      cur_bond->setBrokenBond(intersectWithInitialCrack(cur_bond));
      if (cur_bond->getBrokenBond() == true)
        cur_bond->setPairwiseForce(zero);

    }

    cur_point->calculateLocalDamage();

  }
}



void
PeriCrack::createInitialCrack(std::vector<StateMaterialPointP>& statesArray)
{
  Vector3D zero(0.0, 0.0, 0.0);
  int size = statesArray.size();
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_point = statesArray[i];
    StateBondPArray statesArray = cur_point->getStateBonds();
    for (int j = 0; j < statesArray.size(); j++)
    {
      StateBondP cur_bond = statesArray[j];
      cur_bond->setBrokenBond(intersectWithInitialCrack(cur_bond));
      if (cur_bond->getBrokenBond() == true)
      {
        cur_bond->setInfluenceFunction(0.0);
        cur_bond->setForceVectorState(zero);
      }

    }

    cur_point->calculateStateLocalDamage();

  }
}



bool
PeriCrack::intersectWithInitialCrack(const PeriBondP& bond)
{
  bool intersect = false;
  PeriMaterialPointP center_point = bond->getCenterPoint();
  PeriMaterialPointP second_point = bond->getSecondPoint();   
  bool intersectX = interval(getMinPoint().x(), getMaxPoint().x(), center_point->x(), second_point->x());
  bool intersectY = interval(getMinPoint().y(), getMaxPoint().y(), center_point->y(), second_point->y());
  bool intersectZ = interval(getMinPoint().z(), getMaxPoint().z(), center_point->z(), second_point->z());


  intersect = intersectX && intersectY && intersectZ; 


  return intersect;
}




bool
PeriCrack::intersectWithInitialCrack(const StateBondP& bond)
{
  bool intersect = false;
  StateMaterialPointP center_point = bond->getCenterPoint();
  StateMaterialPointP second_point = bond->getSecondPoint();   
  bool intersectX = interval(getMinPoint().x(), getMaxPoint().x(), center_point->x(), second_point->x());
  bool intersectY = interval(getMinPoint().y(), getMaxPoint().y(), center_point->y(), second_point->y());
  bool intersectZ = interval(getMinPoint().z(), getMaxPoint().z(), center_point->z(), second_point->z());


  intersect = intersectX && intersectY && intersectZ; 


  return intersect;
}





bool
PeriCrack::interval(double Min, double Max, double A, double B)
{
  bool intersection = false;
  if (A == B)
    intersection = ((Min <= B) && (B <= Max));
  else
  {
    double left = (Min - B)/(A - B);
    double right = (Max - B)/(A - B);
    intersection = (((left < 0) && (right < 0)) || ((left > 1) && (right > 1)));
    intersection = !(intersection);
  }
  return intersection; 
 
}

