extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <random>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////// FUNCIONES BL CONTINUA

void clip(vector<double> &sol, int lower, int upper) {
  for (auto &val : sol) {
    if (val < lower) {
      val = lower;
    }
    else if (val > upper) {
      val = upper;
    }
  }
}

void increm_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = 0.2*bias[i]+0.4*(dif[i]+bias[i]);
  }
}

void decrement_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = bias[i]-0.4*(dif[i]+bias[i]);
  }
}

/**
 * Aplica el Solis Wets
 *
 * @param  sol solucion a mejorar.
 * @param fitness fitness de la solución.
 */
template <class Random>
void soliswets(vector<double> &sol, double &fitness, double delta, int maxevals, int lower, int upper, Random &random) {
  const size_t dim = sol.size();
  vector<double> bias (dim), dif (dim), newsol (dim);
  double newfit;
  size_t i;

  int evals = 0;
  int num_success = 0;
  int num_failed = 0;

  while (evals < maxevals) {
    std::uniform_real_distribution<double> distribution(0.0, delta);

    for (i = 0; i < dim; i++) {
      dif[i] = distribution(random);
      newsol[i] = sol[i] + dif[i] + bias[i];
    }

    clip(newsol, lower, upper);
    newfit = cec17_fitness(&newsol[0]);
    evals += 1;

    if (newfit < fitness) {
      sol = newsol;
      fitness = newfit;
      increm_bias(bias, dif);
      num_success += 1;
      num_failed = 0;
    }
    else if (evals < maxevals) {

      for (i = 0; i < dim; i++) {
        newsol[i] = sol[i] - dif[i] - bias[i];
      }

      clip(newsol, lower, upper);
      newfit = cec17_fitness(&newsol[0]);
      evals += 1;

      if (newfit < fitness) {
        sol = newsol;
        fitness = newfit;
        decrement_bias(bias, dif);
        num_success += 1;
        num_failed = 0;
      }
      else {
        for (i = 0; i < dim; i++) {
          bias[i] /= 2;
        }

        num_success = 0;
        num_failed += 1;
      }
    }

    if (num_success >= 5) {
      num_success = 0;
      delta *= 2;
    }
    else if (num_failed >= 3) {
      num_failed = 0;
      delta /= 2;
    }
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

class puebloExplorador{
  private:
    // Vector de double para representar la solución (lugar donde está asentado el pueblo)
    vector<double> asentamiento;

    // Valor para almacenar el valor de la solución
    double valorSol;

    // Valor para representar el grado de confrontación que posee este pueblo
    double gradoConfrontacion;

    // Valor para representar la capacidad de exploración que posee este pueblo
    double gradoExplorador;

    // Valor máximo que de iteraciones que puede pasar el pueblo sin mejora
    int pacienciaPueblo;

    // Valor para almacenar el número de veces que se ha explorado el entorno (sin mejora)
    // Equivale a la idea de que un pueblo realice varias incursiones y no encuentre nada
    int numItersNoMejora = 0;

  public:
    // Constructor de la clase
    // ;
    // dim -> dimensión del problema
    // limitInf -> límite inferior del rango del problema
    // limitSup -> límite superior del rango del problema
    // gen -> pseudogenerador de números aleatorios

    puebloExplorador(int dim, double limitInf, double limitSup, mt19937 &gen){

      // Distribución uniforme para el asentamiento aleatorio
      uniform_real_distribution<> distAsentamiento(limitInf, limitSup);

      // Crear asentamiento en lugar aleatorio
      for (int i = 0; i < dim; i++)
        asentamiento.push_back(distAsentamiento(gen));

      // Inicializar valor de la solución
      valorSol = cec17_fitness(&asentamiento[0]);

      // Distribución uniforme para la generación de características del pueblo
      uniform_real_distribution<> distCaracs(0.0, 1.0);

      // Características aleatorias de los pueblos
      gradoConfrontacion = distCaracs(gen);
      gradoExplorador = distCaracs(gen);

      // Distribución uniforme para la generación de 'pacienciaPueblo'
      // IMPORTANTE: se establece que el valor de la paciencia de cada pueblo será
      // un valor entre 20 * dim y 40 * dim
      uniform_int_distribution<> distPaciencia(dim, 5 * dim);

      // Se genera 'pacienciaPueblo'
      pacienciaPueblo = distPaciencia(gen);
    }
    
    // Getters & setters
    vector<double> getAsentamiento(){
      return asentamiento;
    }

    double getValorSol(){
      return valorSol;
    }

    double getGradoConfrontacion(){
      return gradoConfrontacion;
    }

    double getGradoExplorador(){
      return gradoExplorador;
    }

    int getPacienciaPueblo(){
      return pacienciaPueblo;
    }

    int getNumItersNoMejora(){
      return numItersNoMejora;
    }

    void setAsentamiento(const vector<double> &asentamiento){
      this->asentamiento = asentamiento;
    }

    void setValorSol(double valorSol){
      this->valorSol = valorSol;
    }

    void setGradoConfrontacion(double gradoConfrontacion){
      this->gradoConfrontacion = gradoConfrontacion;
    }

    void setGradoExplorador(double gradoExplorador){
      this->gradoExplorador = gradoExplorador;
    }

    void setPacienciaPueblo(int pacienciaPueblo){
      this->pacienciaPueblo =  pacienciaPueblo;
    }

    void setNumItersNoMejora(int numItersNoMejora){
      this->numItersNoMejora = numItersNoMejora;
    }

    // Métodos útiles de la clase

    // Método que incrementa el número de iteraciones en las que no se ha encontrado mejora en 1 unidad
    void incrementaItersNoMejora(){
      numItersNoMejora++;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Función que calcula la distancia euclídea entre dos pueblos
double distPueblos(int dim, puebloExplorador pueblo1, puebloExplorador pueblo2){
  vector<double> posPueblo1 = pueblo1.getAsentamiento();
  vector<double> posPueblo2 = pueblo2.getAsentamiento();
  
  // Inicializar suma
  double suma = 0;

  for (int i = 0; i < dim; i++){
    double diferencia = posPueblo2[i] - posPueblo1[i];
    suma += diferencia * diferencia;    // s = dif^2 
  }

  return sqrt(suma);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Función que se encarga de unir dos pueblos, todos las características del nuevo pueblo
// se consiguen haciendo la media de las características de los pueblos.
// El pueblo unión es añadido y uno de los dos pueblos se elimina.
// ;
// mundo -> referencia al vector de pueblos
// i y j -> índices de los pueblos a unir
void unirPueblos(vector<puebloExplorador> &mundo, int i, int j){

  // Media de soluciones
  vector<double> solPueblo1 = mundo[i].getAsentamiento();
  vector<double> solPueblo2 = mundo[j].getAsentamiento();
  vector<double> solUnion;

  for (int k = 0; k < solPueblo1.size(); k++)
    solUnion.push_back((solPueblo1[k] + solPueblo2[k]) / 2);

  // Media de grados de confrontación
  double confrontacionUnion = (mundo[i].getGradoConfrontacion() + mundo[j].getGradoConfrontacion()) / 2;

  // Media de grados de exploración
  double exploracionUnion = (mundo[i].getGradoExplorador() + mundo[j].getGradoExplorador()) / 2;

  // Media de paciencias
  int pacienciaUnion = (mundo[i].getPacienciaPueblo() + mundo[j].getPacienciaPueblo()) / 2;

  // Calcular solución del pueblo unión
  double valorSolUnion = cec17_fitness(&solUnion[0]);

  // Todos los cambios se guardan en el pueblo con índice j
  mundo[j].setAsentamiento(solUnion);
  mundo[j].setValorSol(valorSolUnion);
  mundo[j].setGradoConfrontacion(confrontacionUnion);
  mundo[j].setGradoExplorador(exploracionUnion);
  mundo[j].setPacienciaPueblo(pacienciaUnion);
  mundo[j].setNumItersNoMejora(0);  // El número de iteraciones sin mejora se reinicia

  // Se elimina el pueblo de índice i
  mundo.erase(mundo.begin() + i);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Función que se encarga de simular una guerra entre dos pueblos. 
// El pueblo de mejor valor gana. El otro queda eliminado
void guerraPueblos(vector<puebloExplorador> &mundo, int i, int j){
  if(mundo[i].getValorSol() < mundo[j].getValorSol())
    mundo.erase(mundo.begin() + j);
  else
    mundo.erase(mundo.begin() + i);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Función que se encarga de realizar la interactuación entre pueblos.
// Si la suma del grado de confrontación es menor que 0.8 entonces los pueblos se unen.
// En caso contrario el mejor pueblo (con mejor valoración) elimina el peor.
// ;
// mundo -> referencia al vector de pueblos
// i y j -> índices de los pueblos a interactuar
// evalsAxt -> referencia al número de evaluaciones
void interactuarPueblos(vector<puebloExplorador> &mundo, int i, int j, int &evalsAct){

  // Verificar si el grado de confrontación es bajo y entonces se unen los pueblos
  if (mundo[i].getGradoConfrontacion() + mundo[j].getGradoConfrontacion() <= 0.8){
    unirPueblos(mundo, i, j);
    evalsAct++;
  }
  else
    guerraPueblos(mundo, i, j);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Función que crea un vector normalizado de la dimensión indicada.
// Sirve para moverse en el espacio.
vector<double> generaNormal(int dim, double limitInf, double limitSup, mt19937 &gen){
  uniform_real_distribution<> distVector(limitInf, limitSup);

  vector<double> vNormal;
  double suma = 0;

  for(int i = 0; i < dim; i++){
    double valorAleat = distVector(gen);  // Crea valor aleatorio, actualizar suma y añadir a vector
    suma += valorAleat * valorAleat;
    vNormal.push_back(valorAleat);
  }

  double normaV = sqrt(suma); // Cálculo de la norma

  for(int i = 0; i < vNormal.size(); i++) // Normalizar vector
    vNormal[i] = vNormal[i] / normaV;


  return vNormal;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

// Función que multiplica un vector por un escalar
void multiplicaEscalar(vector<double> &v, double escalar){
  for(int i = 0; i < v.size(); i++)
    v[i] *= escalar;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////7

// Función que suma 2 vectores teniendo en cuenta los límtes inferior y superior
vector<double> sumaV(const vector<double> &v1, const vector<double> &v2, double limitInf, double limitSup){
  vector<double> suma = v1;

  for(int i = 0; i < suma.size(); i++){
    if (suma[i] + v2[i] < limitInf)
      suma[i] = limitInf;
    else if (suma[i] + v2[i] > limitSup)
      suma[i] = limitSup;
    else
      suma[i] += v2[i];
  }

  return suma;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////7

// Función encargada de simular la navegación y asentamiento de los pueblos
// dim -> dimensión del problema
// limitInf -> límite inferior del problema
// limitSup -> límite superior del problema
// gen -> pseudogenerador de números aleatorios
void evolucionaPueblos(int dim, double limitInf, double limitSup, mt19937 &gen){

  // En primer lugar, se calcula el número de pueblos
  // IMPORTANTE: se establece que el número de pueblos sea 15 * dim
  int numPueblos = 15 * dim;

  // En segundo lugar, se crean los pueblos
  vector<puebloExplorador> mundo;

  for (int i = 0; i < numPueblos; i++)
    mundo.push_back(puebloExplorador(dim, limitInf, limitSup, gen));
  

  // Se establece la distancia de enfrentamiento (dos pueblos a una distancia menor o igual a esta
  // interactuarán bien con una guerra o bien con una unión)
  // IMPORTANTE: Este parámetro está fijado por mí
  const double distInteractuar = dim * dim / 3.0;

  // Cálculo del número máximo de evaluaciones
  const int maxEvals = 10000 * dim;
  int evalsAct = numPueblos;

  // Llevar la cuenta del número de iteraciones
  int iters = 0;

  // Bucle principal del algoritmo
  while (evalsAct < maxEvals){

    // Cada 10 iteraciones se aplica BL continua 'progresiva', al principio pocas iteraciones y
    // al final muchas
    if(iters % 10 == 0 && iters != 0){

        // Aplicar BL a todos los individuos
        for(int i = 0; i < mundo.size() && evalsAct < maxEvals; i++){
            
            // Cálculo del número de evaluaciones a aplicar
            int evalsAplicar;

            if ((maxEvals - evalsAct) < iters)
                evalsAplicar = maxEvals - evalsAct;
            else
                evalsAplicar = iters;

            // Aplicar BL continua
            vector<double> solPueblo = mundo[i].getAsentamiento();
            double valorSol = mundo[i].getValorSol();

            soliswets(solPueblo, valorSol, 0.2, evalsAplicar, -100.0, 100.0, gen);
            evalsAct += evalsAplicar;

            // Actualizar valores si mejora
            if(valorSol != mundo[i].getValorSol()){
                mundo[i].setAsentamiento(solPueblo);
                mundo[i].setValorSol(valorSol);
                mundo[i].setNumItersNoMejora(0);
            }
        }  // FOR i
    } // if iters

    // En segundo lugar, se mira entre posibles enfrentamientos de pueblos
    for (int i = 0; i < mundo.size() && evalsAct < maxEvals; i++)
      for (int j = i + 1; j < mundo.size() && evalsAct < maxEvals; j++){

        // Si la distancia entre dos pueblos es menor
        if (distPueblos(dim, mundo[i], mundo[j]) <= distInteractuar)
          interactuarPueblos(mundo, i, j, evalsAct);

      } // FIN for j
    
    // En tercer lugar, se crea un vector normal aleatorio en esta iteración
    vector<double> vNormal = generaNormal(dim, limitInf, limitSup, gen);

    // Distribución para la generación de posiciones aleatorias del vector normal
    uniform_int_distribution<> distCambioPos(0, vNormal.size() - 1);

    // Se itera sobre todos los pueblos
    for(int i = 0; i < mundo.size() && evalsAct < maxEvals; i++){
      
      // Ya que generar un vector normal para cada pueblo sería muy costoso, voy
      // intercambiando posiciones aleatorias del vector anteriormmente generado
      int posAleat1 = distCambioPos(gen);
      int posAleat2 = distCambioPos(gen);

      // Intercambio de posiciones en el vector normal
      double cambioVNormal = vNormal[posAleat1];
      vNormal[posAleat1] = vNormal[posAleat2];
      vNormal[posAleat2] = cambioVNormal;

      // Genero un vector distancia a partir del grado de exploración del pueblo
      // y del vector normal.
      // IMPORTANTE: el valor del grado de confrontación se multiplica por 50
      vector<double> vDist = vNormal;
      multiplicaEscalar(vDist, 50 * mundo[i].getGradoExplorador());

      // Genero el nuevo posible asentamiento con la suma del asentamiento actual y la distancia
      // generada
      vector<double> posibleAsentamiento = sumaV(mundo[i].getAsentamiento(), vDist, limitInf, limitSup);

      // Se evalúa el posible asentamiento y se actualiza asentamiento del pueblo
      double nuevoValor = cec17_fitness(&posibleAsentamiento[0]);
      evalsAct++;

      // Si asentamiento es mejor o el número de iteraciones ha alcanzado la paciencia del pueblo, actualizar
      if (nuevoValor < mundo[i].getValorSol()){
        mundo[i].setAsentamiento(posibleAsentamiento);
        mundo[i].setValorSol(nuevoValor);
        mundo[i].setNumItersNoMejora(0);
        mundo[i].setGradoExplorador(mundo[i].getGradoExplorador() / 1.02);  // Por cada mejora del pueblo se va reduciendo el grado explorador
      }
      else if (mundo[i].getNumItersNoMejora() == mundo[i].getPacienciaPueblo() && mundo.size() > 5){
        mundo.erase(mundo.begin() + i);
        i--;
      }
      else
        mundo[i].incrementaItersNoMejora();
    } // FIN for i

    iters++; // Incremento de iteraciones
  } // FIN WHILE
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
  vector<int> dims = {10, 30};
  vector<int> seeds = {3, 45, 68, 1003, 3423, 23943, 95433, 65743, 12344, 86754};

  for(int i = 0; i < dims.size(); i++)
    for(int j = 0; j < seeds.size(); j++){ // Se ejecuta 10 veces, hay 10 semillas

      // Inicializar semilla
      std::mt19937 gen(seeds[j]);

      for(int funcid = 1; funcid <= 30; funcid++){
        cec17_init("mejoraPuebloExpl", funcid, dims[i]);
        evolucionaPueblos(dims[i], -100.0, 100.0, gen);
      }
    }
}
