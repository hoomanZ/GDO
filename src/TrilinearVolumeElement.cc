#include <TrilinearVolumeElement.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <array>


using namespace FiniteElement;

TrilinearVolumeElement::TrilinearVolumeElement()
   : d_node_1(new Node()), d_node_2(new Node()), d_node_3(new Node()),
     d_node_4(new Node()), d_node_5(new Node()), d_node_6(new Node()),
     d_node_7(new Node()), d_node_8(new Node()), noTractionForce(0.0)
{
/*  d_nodes_array{{d_node_1, d_node_2, d_node_3,
                 d_node_4, d_node_5, d_node_6,
                 d_node_7, d_node_8}}; */
  std::fill_n(d_nodes_array.begin(), 1, d_node_1);
  std::fill_n(d_nodes_array.begin()+1, 1, d_node_2);
  std::fill_n(d_nodes_array.begin()+2, 1, d_node_3);
  std::fill_n(d_nodes_array.begin()+3, 1, d_node_4);
  std::fill_n(d_nodes_array.begin()+4, 1, d_node_5);
  std::fill_n(d_nodes_array.begin()+5, 1, d_node_6);
  std::fill_n(d_nodes_array.begin()+6, 1, d_node_7);
  std::fill_n(d_nodes_array.begin()+7, 1, d_node_8);
}

TrilinearVolumeElement::TrilinearVolumeElement(NodeP node1, NodeP node2,
                                               NodeP node3, NodeP node4,
                                               NodeP node5, NodeP node6,
                                               NodeP node7, NodeP node8)
    : d_node_1(node1), d_node_2(node2), d_node_3(node3),
      d_node_4(node4), d_node_5(node5), d_node_6(node6),
      d_node_7(node7), d_node_8(node8), noTractionForce(0.0)
{
/*  d_nodes_array = {d_node_1, d_node_2, d_node_3,
                   d_node_4, d_node_5, d_node_6,
                   d_node_7, d_node_8}; */
  std::fill_n(d_nodes_array.begin(), 1, d_node_1);
  std::fill_n(d_nodes_array.begin()+1, 1, d_node_2);
  std::fill_n(d_nodes_array.begin()+2, 1, d_node_3);
  std::fill_n(d_nodes_array.begin()+3, 1, d_node_4);
  std::fill_n(d_nodes_array.begin()+4, 1, d_node_5);
  std::fill_n(d_nodes_array.begin()+5, 1, d_node_6);
  std::fill_n(d_nodes_array.begin()+6, 1, d_node_7);
  std::fill_n(d_nodes_array.begin()+7, 1, d_node_8);
}

TrilinearVolumeElement::TrilinearVolumeElement(std::array <NodeP, 8> nodes)
    : d_node_1(nodes[0]), d_node_2(nodes[1]), d_node_3(nodes[2]),
      d_node_4(nodes[3]), d_node_5(nodes[4]), d_node_6(nodes[5]),
      d_node_7(nodes[6]), d_node_8(nodes[7]), noTractionForce(0.0)
{
/*  d_nodes_array = {d_node_1, d_node_2, d_node_3,
                   d_node_4, d_node_5, d_node_6,
                   d_node_7, d_node_8}; */
  std::fill_n(d_nodes_array.begin(), 1, d_node_1);
  std::fill_n(d_nodes_array.begin()+1, 1, d_node_2);
  std::fill_n(d_nodes_array.begin()+2, 1, d_node_3);
  std::fill_n(d_nodes_array.begin()+3, 1, d_node_4);
  std::fill_n(d_nodes_array.begin()+4, 1, d_node_5);
  std::fill_n(d_nodes_array.begin()+5, 1, d_node_6);
  std::fill_n(d_nodes_array.begin()+6, 1, d_node_7);
  std::fill_n(d_nodes_array.begin()+7, 1, d_node_8);
}

TrilinearVolumeElement::TrilinearVolumeElement(std::vector <NodeP> nodes,
                                               int n1, int n2, int n3, int n4,
                                               int n5, int n6, int n7, int n8)
     :d_node_1(nodeFunction(nodes, n1)), d_node_2(nodeFunction(nodes, n2)),
      d_node_3(nodeFunction(nodes, n3)), d_node_4(nodeFunction(nodes, n4)),
      d_node_5(nodeFunction(nodes, n5)), d_node_6(nodeFunction(nodes, n6)),
      d_node_7(nodeFunction(nodes, n7)), d_node_8(nodeFunction(nodes, n8)), 
      noTractionForce(0.0)
{
  std::fill_n(d_nodes_array.begin(), 1, d_node_1);
  std::fill_n(d_nodes_array.begin()+1, 1, d_node_2);
  std::fill_n(d_nodes_array.begin()+2, 1, d_node_3);
  std::fill_n(d_nodes_array.begin()+3, 1, d_node_4);
  std::fill_n(d_nodes_array.begin()+4, 1, d_node_5);
  std::fill_n(d_nodes_array.begin()+5, 1, d_node_6);
  std::fill_n(d_nodes_array.begin()+6, 1, d_node_7);
  std::fill_n(d_nodes_array.begin()+7, 1, d_node_8);  
}


TrilinearVolumeElement::~TrilinearVolumeElement()
{
}

void
TrilinearVolumeElement::setArray(NodeP node1, NodeP node2,
                                 NodeP node3, NodeP node4,
                                 NodeP node5, NodeP node6,
                                 NodeP node7, NodeP node8)

{
    d_node_1 = node1;
    d_node_2 = node2;
    d_node_3 = node3;
    d_node_4 = node4;
    d_node_5 = node5;
    d_node_6 = node6;
    d_node_7 = node7;
    d_node_8 = node8;
/*    d_nodes_array = {d_node_1, d_node_2, d_node_3,
                     d_node_4, d_node_5, d_node_6,
                     d_node_7, d_node_8}; */
  std::fill_n(d_nodes_array.begin(), 1, d_node_1);
  std::fill_n(d_nodes_array.begin()+1, 1, d_node_2);
  std::fill_n(d_nodes_array.begin()+2, 1, d_node_3);
  std::fill_n(d_nodes_array.begin()+3, 1, d_node_4);
  std::fill_n(d_nodes_array.begin()+4, 1, d_node_5);
  std::fill_n(d_nodes_array.begin()+5, 1, d_node_6);
  std::fill_n(d_nodes_array.begin()+6, 1, d_node_7);
  std::fill_n(d_nodes_array.begin()+7, 1, d_node_8);
  noTractionForce.set(0.0);
}

void
TrilinearVolumeElement::setArray(const std::array <NodeP, 8>& nodesArray)

{
    d_node_1 = nodesArray[0];
    d_node_2 = nodesArray[1];
    d_node_3 = nodesArray[2];
    d_node_4 = nodesArray[3];
    d_node_5 = nodesArray[4];
    d_node_6 = nodesArray[5];
    d_node_7 = nodesArray[6];
    d_node_8 = nodesArray[7];
/*    d_nodes_array = {d_node_1, d_node_2, d_node_3,
                     d_node_4, d_node_5, d_node_6,
                     d_node_7, d_node_8}; */
  std::fill_n(d_nodes_array.begin(), 1, d_node_1);
  std::fill_n(d_nodes_array.begin()+1, 1, d_node_2);
  std::fill_n(d_nodes_array.begin()+2, 1, d_node_3);
  std::fill_n(d_nodes_array.begin()+3, 1, d_node_4);
  std::fill_n(d_nodes_array.begin()+4, 1, d_node_5);
  std::fill_n(d_nodes_array.begin()+5, 1, d_node_6);
  std::fill_n(d_nodes_array.begin()+6, 1, d_node_7);
  std::fill_n(d_nodes_array.begin()+7, 1, d_node_8);
  noTractionForce.set(0.0);
}

void 
TrilinearVolumeElement::operator = (const TrilinearVolumeElement& triElement)
  
{
  setNodeP1(triElement.getNodeP1());
  setNodeP2(triElement.getNodeP2());
  setNodeP3(triElement.getNodeP3());
  setNodeP4(triElement.getNodeP4());
  setNodeP5(triElement.getNodeP5());
  setNodeP6(triElement.getNodeP6());
  setNodeP7(triElement.getNodeP7());
  setNodeP8(triElement.getNodeP8());
  setArray(triElement.getArray());
  noTractionForce.set(0.0);
}

void 
TrilinearVolumeElement::operator = (const TrilinearVolumeElementSP& triElement)
{
  setNodeP1(triElement->getNodeP1());
  setNodeP2(triElement->getNodeP2());
  setNodeP3(triElement->getNodeP3());
  setNodeP4(triElement->getNodeP4());
  setNodeP5(triElement->getNodeP5());
  setNodeP6(triElement->getNodeP6());
  setNodeP7(triElement->getNodeP7());
  setNodeP8(triElement->getNodeP8());
  setArray(triElement->getArray());
//  noTractionForce.set(0.0);
}


NodeP
TrilinearVolumeElement::nodeFunction(std::vector <NodeP> nodes, int i)
{
  bool isNumberAnID = false;
  for (auto iter = nodes.begin(); iter != nodes.end(); iter++)
  {
    NodeP node = *iter;
    if (i == node->getID())
    {
      return node;
      isNumberAnID = true;
      break;
    }
  }
  if (!isNumberAnID)
  {
    std::cout << " ****WARNING***** The ID is invalid. A node is defined by default." << std::endl;
    NodeP defaultNode(new Node());
    return defaultNode;
  }
}


//                             interval=[0, 1]
double
TrilinearVolumeElement::testFunction(int a, double xi1, 
                                     double xi2, double xi3)
{
   switch (a)
   {
     case 1:
     {
       return (1-xi1)*(1-xi2)*(1-xi3);
       break;
     }
     case 2:
     {
       return xi1*(1-xi2)*(1-xi3);
       break;
     }
     case 3:
     {
       return xi1*xi2*(1-xi3);
       break;
     }
     case 4:
     {
       return (1-xi1)*xi2*(1-xi3);
       break;
     }
     case 5:
     {
       return (1-xi1)*xi2*xi3;
       break;
     }
     case 6:
     {
       return (1-xi1)*(1-xi2)*xi3;
       break;
     }
     case 7:
     {
       return xi1*(1-xi2)*xi3;
       break;
     }
     case 8:
     {
       return xi1*xi2*xi3;
       break;
     }
     default:
     {
       std::cout << " The index a shows the number of local nodes in the element which is from 1 to 8.";
       break;
     }
  } //end of switch
}//end of function


double
TrilinearVolumeElement::derivativeTestFunctionToXi(int a, int j,
                                       double xi1, double xi2, double xi3)
{
  switch (a)
  {
     case 1:
     {
       switch (j)
       {
          case 1:
          {
             return -(1-xi2)*(1-xi3);
             break;
          }       
          case 2:
          {
             return -(1-xi1)*(1-xi3);
             break;
          }       
          case 3:
          {
             return -(1-xi1)*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }      
     case 2:
     {
       switch (j)
       {
          case 1:
          {
             return (1-xi2)*(1-xi3);
             break;
          }       
          case 2:
          {
             return -xi1*(1-xi3);
             break;
          }       
          case 3:
          {
             return -xi1*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 3:
     {
       switch (j)
       {
          case 1:
          {
             return xi2*(1-xi3);
             break;
          }       
          case 2:
          {
             return xi1*(1-xi3);
             break;
          }       
          case 3:
          {
             return -xi1*xi2;
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 4:
     {
       switch (j)
       {
          case 1:
          {
             return -xi2*(1-xi3);
             break;
          }       
          case 2:
          {
             return (1-xi1)*(1-xi3);
             break;
          }       
          case 3:
          {
             return -(1-xi1)*xi2;
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 5:
     {
       switch (j)
       {
          case 1:
          {
             return -xi2*xi3;
             break;
          }       
          case 2:
          {
             return (1-xi1)*xi3;
             break;
          }       
          case 3:
          {
             return (1-xi1)*xi2;
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 6:
     {
       switch (j)
       {
          case 1:
          {
             return -(1-xi2)*xi3;
             break;
          }       
          case 2:
          {
             return -(1-xi1)*xi3;
             break;
          }       
          case 3:
          {
             return (1-xi1)*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 7:
     {
       switch (j)
       {
          case 1:
          {
             return (1-xi2)*xi3;
             break;
          }       
          case 2:
          {
             return -xi1*xi3;
             break;
          }       
          case 3:
          {
             return xi1*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 8:
     {
       switch (j)
       {
          case 1:
          {
             return xi2*xi3;
             break;
          }       
          case 2:
          {
             return xi1*xi3;
             break;
          }       
          case 3:
          {
             return xi1*xi2;
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }
     default:
     {
       std::cout << " The index a shows the number of local nodes in the element which is from 1 to 8.";
       break;
     }
   }//end of switch (a) 
}//end of function


Matrix3D
TrilinearVolumeElement::jacobianMatrix(double xi1, double xi2, double xi3) //derivativePositionToXi             
{
  double sum = 0.0;
  Matrix3D jacob(0.0);
  for(int i = 1; i < 4; i++) {
    for(int j = 1; j < 4; j++) {
      sum = 0.0;
      for (int a = 1; a < 9; a++) {
        sum = sum + (d_nodes_array[a-1])->operator()(i)*derivativeTestFunctionToXi(a,j,xi1,xi2,xi3);
        /*if (i != j) {
        std::cout << "i= " << i << ",  j= " << j << ",  a= " << a << std::endl;
        std::cout << "d_nodes_array[" << a-1 << "]= " << (d_nodes_array[a-1])->operator()(i) << std::endl;
        std::cout << "derivativeTestFunctionToXi(" << a << "," << j << "," << xi1 << ","
                  << xi2 << "," << xi3 << ")= " << derivativeTestFunctionToXi(a,j,xi1,xi2,xi3) << std::endl;
        }*/
//        sum = sum + (d_nodes_array[a-1])->(i)*derivativeTestFunctionToXi(a,j,xi1,xi2,xi3);
      }
      jacob(i-1, j-1)= sum;
    }
  }
  return jacob;
}  


double
TrilinearVolumeElement::derivativeTestFunctionToPosition(int a, int i,
                                                         double xi1, double xi2, double xi3)
{
  Matrix3D invJacob = jacobianMatrix(xi1, xi2, xi3).Inverse();
  double sum = 0.0;
  for(int j = 1; j < 4; j++) {
    sum = sum + derivativeTestFunctionToXi(a,j,xi1,xi2,xi3)*invJacob(j-1, i-1);
  }
  return sum;
}

  
//                           interval= [-1, 1]

double
TrilinearVolumeElement::symTestFunction(int a, double xi1, 
                                        double xi2, double xi3)
{
   switch (a)
   {
     case 1:
     {
       return 0.125*(1-xi1)*(1-xi2)*(1-xi3);
       break;
     }
     case 2:
     {
       return 0.125*(1+xi1)*(1-xi2)*(1-xi3);
       break;
     }
     case 3:
     {
       return 0.125*(1+xi1)*(1+xi2)*(1-xi3);
       break;
     }
     case 4:
     {
       return 0.125*(1-xi1)*(1+xi2)*(1-xi3);
       break;
     }
     case 5:
     {
       return 0.125*(1-xi1)*(1+xi2)*(1+xi3);
       break;
     }
     case 6:
     {
       return 0.125*(1-xi1)*(1-xi2)*(1+xi3);
       break;
     }
     case 7:
     {
       return 0.125*(1+xi1)*(1-xi2)*(1+xi3);
       break;
     }
     case 8:
     {
       return 0.125*(1+xi1)*(1+xi2)*(1+xi3);
       break;
     }
     default:
     {
       std::cout << " The index a shows the number of local nodes in the element which is from 1 to 8.";
       break;
     }
  } //end of switch
}//end of function

double
TrilinearVolumeElement::symDerivativeTestFunctionToXi(int a, int j,
                                       double xi1, double xi2, double xi3)
{
  switch (a)
  {
     case 1:
     {
       switch (j)
       {
          case 1:
          {
             return -0.125*(1-xi2)*(1-xi3);
             break;
          }       
          case 2:
          {
             return -0.125*(1-xi1)*(1-xi3);
             break;
          }       
          case 3:
          {
             return -0.125*(1-xi1)*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }      
     case 2:
     {
       switch (j)
       {
          case 1:
          {
             return 0.125*(1-xi2)*(1-xi3);
             break;
          }       
          case 2:
          {
             return -0.125*(1+xi1)*(1-xi3);
             break;
          }       
          case 3:
          {
             return -0.125*(1+xi1)*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 3:
     {
       switch (j)
       {
          case 1:
          {
             return 0.125*(1+xi2)*(1-xi3);
             break;
          }       
          case 2:
          {
             return 0.125*(1+xi1)*(1-xi3);
             break;
          }       
          case 3:
          {
             return -0.125*(1+xi1)*(1+xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 4:
     {
       switch (j)
       {
          case 1:
          {
             return -0.125*(1+xi2)*(1-xi3);
             break;
          }       
          case 2:
          {
             return 0.125*(1-xi1)*(1-xi3);
             break;
          }       
          case 3:
          {
             return -0.125*(1-xi1)*(1+xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 5:
     {
       switch (j)
       {
          case 1:
          {
             return -0.125*(1+xi2)*(1+xi3);
             break;
          }       
          case 2:
          {
             return 0.125*(1-xi1)*(1+xi3);
             break;
          }       
          case 3:
          {
             return 0.125*(1-xi1)*(1+xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 6:
     {
       switch (j)
       {
          case 1:
          {
             return -0.125*(1-xi2)*(1+xi3);
             break;
          }       
          case 2:
          {
             return -0.125*(1-xi1)*(1+xi3);
             break;
          }       
          case 3:
          {
             return 0.125*(1-xi1)*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 7:
     {
       switch (j)
       {
          case 1:
          {
             return 0.125*(1-xi2)*(1+xi3);
             break;
          }       
          case 2:
          {
             return -0.125*(1+xi1)*(1+xi3);
             break;
          }       
          case 3:
          {
             return 0.125*(1+xi1)*(1-xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }              
     case 8:
     {
       switch (j)
       {
          case 1:
          {
             return 0.125*(1+xi2)*(1+xi3);
             break;
          }       
          case 2:
          {
             return 0.125*(1+xi1)*(1+xi3);
             break;
          }       
          case 3:
          {
             return 0.125*(1+xi1)*(1+xi2);
             break;
          }
          default:
          {
             std::cout << " The index j indicates coordinates x, y, and z so j is from 1 to 3";
             break;
          }
       };
       break;
     }
     default:
     {
       std::cout << " The index a shows the number of local nodes in the element which is from 1 to 8.";
       break;
     }
   }//end of switch (a) 
}//end of function


Matrix3D
TrilinearVolumeElement::symJacobianMatrix(double xi1, double xi2, double xi3) //derivativePositionToXi             
{
  double sum = 0.0;
  Matrix3D jacob(0.0);
  for(int i = 1; i < 4; i++) {
    for(int j = 1; j < 4; j++) {
      sum = 0.0;
      for (int a = 1; a < 9; a++) {
        sum = sum + (d_nodes_array[a-1])->operator()(i)*symDerivativeTestFunctionToXi(a,j,xi1,xi2,xi3);
        /*if (i != j) {
        std::cout << "i= " << i << ",  j= " << j << ",  a= " << a << std::endl;
        std::cout << "d_nodes_array[" << a-1 << "]= " << (d_nodes_array[a-1])->operator()(i) << std::endl;
        std::cout << "symDerivativeTestFunctionToXi(" << a << "," << j << "," << xi1 << ","
                  << xi2 << "," << xi3 << ")= " << symDerivativeTestFunctionToXi(a,j,xi1,xi2,xi3) << std::endl;
        }*/
//        sum = sum + (d_nodes_array[a-1])->(i)*symDerivativeTestFunctionToXi(a,j,xi1,xi2,xi3);
      }
      jacob(i-1, j-1)= sum;
    }
  }
  return jacob;
}  


double
TrilinearVolumeElement::symDerivativeTestFunctionToPosition(int a, int i,
                                                         double xi1, double xi2, double xi3)
{
  Matrix3D invJacob = symJacobianMatrix(xi1, xi2, xi3).Inverse();
  double sum = 0.0;
  for(int j = 1; j < 4; j++) {
    sum = sum + symDerivativeTestFunctionToXi(a,j,xi1,xi2,xi3)*invJacob(j-1, i-1);
  }
  return sum;
}
 


double
TrilinearVolumeElement::symXtoXi1(const double& x)
{
  double numerator = x-d_node_1->x();
  double denominator = d_node_2->x()-d_node_1->x();
  return (-1+2*numerator/denominator);
}
double
TrilinearVolumeElement::symYtoXi2(const double& y)
{
  double numerator = y-d_node_1->y();
  double denominator = d_node_4->y()-d_node_1->y();
  return (-1+2*numerator/denominator);
}
double
TrilinearVolumeElement::symZtoXi3(const double& z)
{
  double numerator = z-d_node_1->z();
  double denominator = d_node_6->z()-d_node_1->z();
  return (-1+2*numerator/denominator);
}



double
TrilinearVolumeElement::symXi1toX(const double& xi1)
{
  double denominator = d_node_2->x()-d_node_1->x();
  return (d_node_1->x()+0.5*denominator*(1+xi1));
}
double
TrilinearVolumeElement::symXi2toY(const double& xi2)
{
  double denominator = d_node_4->y()-d_node_1->y();
  return (d_node_1->y()+0.5*denominator*(1+xi2));
}
double
TrilinearVolumeElement::symXi3toZ(const double& xi3)
{
  double denominator = d_node_6->z()-d_node_1->z();
  return (d_node_1->z()+0.5*denominator*(1+xi3));
}




double
TrilinearVolumeElement::symInterpolatedValue(const std::array<double, 8>& nodesValues, 
                                             const double& xi1, const double& xi2, const double& xi3)
{
  double sum = 0.0;
  for(int a = 1; a < 9; a++) {
    sum += nodesValues[a-1]*symTestFunction(a, xi1, xi2, xi3);
  }
  return sum;
}


Vector3D
TrilinearVolumeElement::symInterpolatedVector(const std::array<Vector3D, 8>& nodesVectors,
                                  const double& xi1, const double& xi2, const double& xi3)
{
  Vector3D interpolatedVector(0, 0, 0);
  std::array<double, 8> arrayXvec;
  std::array<double, 8> arrayYvec;
  std::array<double, 8> arrayZvec;

  std::fill_n(arrayXvec.begin(), 1, nodesVectors[0].x());
  std::fill_n(arrayXvec.begin()+1, 1, nodesVectors[1].x());
  std::fill_n(arrayXvec.begin()+2, 1, nodesVectors[2].x());
  std::fill_n(arrayXvec.begin()+3, 1, nodesVectors[3].x());
  std::fill_n(arrayXvec.begin()+4, 1, nodesVectors[4].x());
  std::fill_n(arrayXvec.begin()+5, 1, nodesVectors[5].x());
  std::fill_n(arrayXvec.begin()+6, 1, nodesVectors[6].x());
  std::fill_n(arrayXvec.begin()+7, 1, nodesVectors[7].x()); 

  std::fill_n(arrayYvec.begin(), 1, nodesVectors[0].y());
  std::fill_n(arrayYvec.begin()+1, 1, nodesVectors[1].y());
  std::fill_n(arrayYvec.begin()+2, 1, nodesVectors[2].y());
  std::fill_n(arrayYvec.begin()+3, 1, nodesVectors[3].y());
  std::fill_n(arrayYvec.begin()+4, 1, nodesVectors[4].y());
  std::fill_n(arrayYvec.begin()+5, 1, nodesVectors[5].y());
  std::fill_n(arrayYvec.begin()+6, 1, nodesVectors[6].y());
  std::fill_n(arrayYvec.begin()+7, 1, nodesVectors[7].y());

  std::fill_n(arrayZvec.begin(), 1, nodesVectors[0].z());
  std::fill_n(arrayZvec.begin()+1, 1, nodesVectors[1].z());
  std::fill_n(arrayZvec.begin()+2, 1, nodesVectors[2].z());
  std::fill_n(arrayZvec.begin()+3, 1, nodesVectors[3].z());
  std::fill_n(arrayZvec.begin()+4, 1, nodesVectors[4].z());
  std::fill_n(arrayZvec.begin()+5, 1, nodesVectors[5].z());
  std::fill_n(arrayZvec.begin()+6, 1, nodesVectors[6].z());
  std::fill_n(arrayZvec.begin()+7, 1, nodesVectors[7].z());

  interpolatedVector.x(symInterpolatedValue(arrayXvec, xi1, xi2, xi3)); 
  interpolatedVector.y(symInterpolatedValue(arrayYvec, xi1, xi2, xi3));
  interpolatedVector.z(symInterpolatedValue(arrayZvec, xi1, xi2, xi3));
   
  return interpolatedVector;
}


Vector3D
TrilinearVolumeElement::symNormVecXi3Max(const double& xi1, const double& xi2)
{
  Matrix3D JM = symJacobianMatrix(xi1, xi2, 1);
  Vector3D vec0(JM(0, 0), JM(1, 0), JM(2, 0));
  Vector3D vec1(JM(0, 1), JM(1, 1), JM(2, 1));
  return vec0.cross(vec1);
}
Vector3D
TrilinearVolumeElement::symNormVecXi3Min(const double& xi1, const double& xi2)
{
  Matrix3D JM = symJacobianMatrix(xi1, xi2, -1);
  Vector3D vec0(JM(0, 0), JM(1, 0), JM(2, 0));
  Vector3D vec1(JM(0, 1), JM(1, 1), JM(2, 1));
  return vec1.cross(vec0);
}
Vector3D
TrilinearVolumeElement::symNormVecXi2Max(const double& xi1, const double& xi3)
{
  Matrix3D JM = symJacobianMatrix(xi1, 1, xi3);
  Vector3D vec0(JM(0, 0), JM(1, 0), JM(2, 0));
  Vector3D vec1(JM(0, 2), JM(1, 2), JM(2, 2));
  return vec1.cross(vec0);
}
Vector3D
TrilinearVolumeElement::symNormVecXi2Min(const double& xi1, const double& xi3)
{
  Matrix3D JM = symJacobianMatrix(xi1, -1, xi3);
  Vector3D vec0(JM(0, 0), JM(1, 0), JM(2, 0));
  Vector3D vec1(JM(0, 2), JM(1, 2), JM(2, 2));
  return vec0.cross(vec1);
}
Vector3D
TrilinearVolumeElement::symNormVecXi1Max(const double& xi2, const double& xi3)
{
  Matrix3D JM = symJacobianMatrix(1, xi2, xi3);
  Vector3D vec0(JM(0, 1), JM(1, 1), JM(2, 1));
  Vector3D vec1(JM(0, 2), JM(1, 2), JM(2, 2));
  return vec0.cross(vec1);
}
Vector3D
TrilinearVolumeElement::symNormVecXi1Min(const double& xi2, const double& xi3)
{
  Matrix3D JM = symJacobianMatrix(-1, xi2, xi3);
  Vector3D vec0(JM(0, 1), JM(1, 1), JM(2, 1));
  Vector3D vec1(JM(0, 2), JM(1, 2), JM(2, 2));
  return vec1.cross(vec0);
}


std::array<double, 8>
TrilinearVolumeElement::densityArray() {
  std::array<double, 8> arrayDensity;
  std::fill_n(arrayDensity.begin(), 1,   d_node_1->densityNode());
  std::fill_n(arrayDensity.begin()+1, 1, d_node_2->densityNode());
  std::fill_n(arrayDensity.begin()+2, 1, d_node_3->densityNode());
  std::fill_n(arrayDensity.begin()+3, 1, d_node_4->densityNode());
  std::fill_n(arrayDensity.begin()+4, 1, d_node_5->densityNode());
  std::fill_n(arrayDensity.begin()+5, 1, d_node_6->densityNode());
  std::fill_n(arrayDensity.begin()+6, 1, d_node_7->densityNode());
  std::fill_n(arrayDensity.begin()+7, 1, d_node_8->densityNode());   
  return arrayDensity;
}


bool
TrilinearVolumeElement::face(const int& i, const int& j,
                             const int& k, const int& l)
{
  noTractionForce.set(0);
  if (!(d_nodes_array[i-1]->surfaceForce() == noTractionForce) && 
      !(d_nodes_array[j-1]->surfaceForce() == noTractionForce) && 
      !(d_nodes_array[k-1]->surfaceForce() == noTractionForce) && 
      !(d_nodes_array[l-1]->surfaceForce() == noTractionForce) )
      return true;
  else return false;
}


std::array<Vector3D, 8>
TrilinearVolumeElement::surfaceForceArray(const int& i, const int& j,
                                          const int& k, const int& l)
{
  noTractionForce.set(0.0);
  std::array<Vector3D, 8> surfaceArray;
  for (int index = 0; index < 8; index++) {
    if ((index != i-1) && (index != j-1) && (index != k-1) && (index != l-1)) {
        std::fill_n(surfaceArray.begin()+index, 1, noTractionForce);
    }
  }  
  std::fill_n(surfaceArray.begin()+i-1, 1, d_nodes_array[i-1]->surfaceForce());
  std::fill_n(surfaceArray.begin()+j-1, 1, d_nodes_array[j-1]->surfaceForce());
  std::fill_n(surfaceArray.begin()+k-1, 1, d_nodes_array[k-1]->surfaceForce());
  std::fill_n(surfaceArray.begin()+l-1, 1, d_nodes_array[l-1]->surfaceForce());
  return surfaceArray;
}



std::array<Vector3D, 8>
TrilinearVolumeElement::bodyForceArray()                                          
{
  std::array<Vector3D, 8> bodyArray;
  for (int index = 0; index < 8; index++) {
     std::fill_n(bodyArray.begin()+index, 1, d_nodes_array[index]->bodyForce());
  }  
  return bodyArray;
}


std::array<Vector3D, 8>
TrilinearVolumeElement::oldDisplacementArray()                                          
{
  std::array<Vector3D, 8> oldDispArray;
  for (int index = 0; index < 8; index++) {
     std::fill_n(oldDispArray.begin()+index, 1, d_nodes_array[index]->oldDisplacement());
  }  
  return oldDispArray;
}


std::array<Vector3D, 8>
TrilinearVolumeElement::displacementArray()                                          
{
  std::array<Vector3D, 8> dispArray;
  for (int index = 0; index < 8; index++) {
     std::fill_n(dispArray.begin()+index, 1, d_nodes_array[index]->displacement());
  }  
  return dispArray;
}


std::array<Vector3D, 8>
TrilinearVolumeElement::newDisplacementArray()                                          
{
  std::array<Vector3D, 8> newDispArray;
  for (int index = 0; index < 8; index++) {
     std::fill_n(newDispArray.begin()+index, 1, d_nodes_array[index]->newDisplacement());
  }  
  return newDispArray;
}


double
TrilinearVolumeElement::twoDJacobian(const int& i, const int& j,
                                     const int& k, const int& l,
                                     const double& xi1, const double& xi2)
{
  if ((i == 5) && (j == 6) && (k == 7) && (l == 8)) {
    return symNormVecXi3Max(xi1, xi2).lengthSq();
  } else if ((i == 1) && (j == 2 ) && (k == 3) && (l == 4)) {
    return symNormVecXi3Min(xi1, xi2).lengthSq();
  } else if ((i == 3) && (j == 4) && (k == 5) && (l == 8)) {
    return symNormVecXi2Max(xi1, xi2).lengthSq();
  } else if ((i == 1) && (j == 2) && (k == 6) && (l == 7)) {
    return symNormVecXi2Min(xi1, xi2).lengthSq();
  } else if ((i == 2) && (j == 3) && (k == 7) && (l == 8)) {
    return symNormVecXi1Max(xi1, xi2).lengthSq();
  } else if ((i == 1) && (j == 4) && (k == 5) && (l == 6)) {
    return symNormVecXi1Min(xi1, xi2).lengthSq();
  } else {
           std::cout << " The possible amounts of i, j, k, and l are:" << std::endl;
           std::cout << " (i, j, k, l) = (5, 6, 7, 8) for plane xi3 = 1" << std::endl;
           std::cout << " (i, j, k, l) = (1, 2, 3, 4) for plane xi3 = -1" << std::endl;
           std::cout << " (i, j, k, l) = (3, 4, 5, 8) for plane xi2 = 1" << std::endl;
           std::cout << " (i, j, k, l) = (1, 2, 6, 7) for plane xi2 = -1" << std::endl;
           std::cout << " (i, j, k, l) = (2, 3, 7, 8) for plane xi1 = 1" << std::endl;
           std::cout << " (i, j, k, l) = (1, 4, 5, 6) for plane xi1 = -1" << std::endl;
           return -1;
  }
}                            




                         
