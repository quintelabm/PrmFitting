#include "HCV_Model.h"

using namespace std;

int main(){
  HCV_Model* model = new HCV_Model();
  
  model->solve();//solve
  system("/bin/python3 $(pwd)/../plotDados_e_Modelo.py");
  return 0;
}
