#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
using namespace std;

class Vector {
    int dim;
    double *V;

public:
    friend class Matrix;

    Vector (int d){
        this->dim=d;
        this->V= new double [d];
        for (int i=0; i<d ; i++)
        {
            this->V[i]=0;
        }

    }


    Vector (int d, double *X)
    {
        this->dim= d;
        this->V= new double [d];
        for(int i=0; i<d; i++)
        {
            this->V[i]=X[i];
        }
    }

    Vector (const Vector &X)
    {
        this->dim=X.dim;
        this->V= new double [this->dim];
        for (int i=0; i<this->dim ; i++)
        {
            this->V[i]=X.V[i];
        }
    }


    ~ Vector (){
        delete[] this->V ;
        this->dim=0;
    }


    Vector& operator= (const Vector& R){

        this->dim = R.dim;
        V= new double [this->dim];
        for(int i=0; i<(this->dim); i++){
            this->V[i] = R.V[i];
        }
        return *this;
    }

    bool operator== (const Vector& R){
        int T=1;
        for(int i=0; i<(this->dim); i++){

            if (this->V[i] != R.V[i]){
                T=0;
            }
        }

        if(T==0){
            return false;
        }else{
            return true;
        }
    }

    bool operator!= (const Vector& R){
        int T=0;
        for(int i=0; i<(this->dim); i++){
            if (this->V[i] != R.V[i]){
                T=1;
            }
        }

        if(T==0){
            return false;
        }else{
            return true;
        }
    }

    Vector operator+ (Vector R){
        Vector tmp (R);

        for(int i=0; i<R.dim; i++)
        {
            tmp.V[i]+= this->V[i];
        }
        return tmp;
    }

    Vector abs (){
        for (int i=0; i<this->dim ; i++)
        {
            this->V[i]= fabs(this->V[i]);
        }

        return *this;
    }

    void print (){
        cout << endl << "Vector = ";
        for(int i=0; i< this->dim; i++) { cout<< "x" << i+1 << " = " <<this->V[i]<<"\n"; }
        cout << endl;
    }

    friend Vector operator- (Vector L, Vector R){
        Vector tmp (L);

        for(int i = 0; i < tmp.dim; i++) { tmp.V[i]-= R.V[i]; }

        return tmp;
    }
};

class Matrix{
    int dim;
    double *M;

public:

    Matrix (int d){
        this->dim=d;
        this->M= new double [d*d];
        for(int i=0; i<(d*d) ; i++){
            this->M[i]=0;
        }
    }

    Matrix (int d, double* X){
        this->dim=d;
        this->M= new double [d*d];
        for(int i=0; i< d * d ; i++){
            this->M[i] = X[i];
        }
    }

    Matrix (const Matrix& r){
        this->dim = r.dim;
        this->M = new double [dim*dim];
        for(int i=0; i<(dim*dim) ; i++){
            this->M[i] = r.M[i];
        }
    }

    ~ Matrix (){
        delete[] this->M;
        this->dim=0;
    }


    Matrix& operator= (const Matrix& R){

        this->dim = R.dim;

        this->M = new double [dim*dim];

        for(int i=0; i<(dim*dim); i++){
            this->M[i] = R.M[i];
        }

        return *this;
    }

    Matrix operator* (const Matrix &R){

        Matrix tmp(R.dim);

        for (int i = 0; i < this->dim; i++){

            for (int j = 0; j < this->dim; j++){

                for (int k = 0; k < this->dim; k++){


                    tmp.M[(this->dim * i) + j] += (this->M[(this->dim * i) + k] * R.M[(this->dim * k) + j]);

                }
            }
        }

        return tmp;

    }

    Vector operator* (Vector &R){

        Vector X (R.dim);

        for(int l=0; l<dim; l++){

            for(int i=(l*dim) , j=0 ; i<((l+1)*dim) ; j++, i++){

                X.V[l] += ( (this->M[i]) * R.V[j] ) ;

            }
        }

        return X;
    }


    void print (){

        for (int i = 0; i < this->dim; i++){

            for (int j = 0; j < this->dim; j++){

                cout << this->M[(this->dim * i) + j] << " ";

            }

            cout << endl;
        }

        cout << endl;
    }

    friend Matrix operator- (Matrix L, Matrix R){

        Matrix tmp (L);

        for(int i=0; i< R.dim * R.dim ; i++) {

            tmp.M[i]-= R.M[i];

        }

        return tmp;
    }

};


int main(){
    cout<< setprecision(6)<< fixed;
    srand(time(NULL));

    int dim;
    cout<<"Enter the dim of arrays and Matrices "<<endl;
    cin>>dim;

    double M1 [dim*dim];
    for(int i = 0; i < dim*dim; i++){

        M1[i] = (rand() % 10000)/100.0 -50;

    }

    for(int i = 0; i < dim; i++){

        M1[i*(dim+1)] *= 100;
    }
    Matrix A (dim, M1);
    A.print();
    cout<<endl;

    double V1 [dim];

    for(int i = 0; i < dim; i++){
        V1[i] = ((i+1) * 10) - 1;
    }

    Vector b (dim, V1);
    b.print();
    cout<<endl;

    double M2 [dim * dim];
    for(int i=0; i<dim * dim; i++){

        if((i % ( dim + 1 )) == 0){

            M2[i] = 1 / M1[i];

        }

        else{

            M2[i]=0;

        }

    }
    Matrix D (dim, M2);
    D.print();
    cout<<endl;

    double M3 [dim * dim];
    for(int i = 0; i < dim * dim; i++){
        if((i % (dim + 1)) == 0){

            M3[i] = 1;

        }

        else{

            M3[i] = 0;
        }

    }
    Matrix I (dim, M3);
    I.print();
    cout<<endl;

    Vector x (dim);
    Vector y (dim);


    for( int i = 0; i < 1000; i++){
        cout << "Iteration [" << i+1 << "]" << endl;
        y = (((I - (A*D)) * x) + (D*b));

        cout << "X = ";
        x.print();

        cout<<endl << "A * x = ";
        (A*x).print();
        cout << endl;

        cout << endl << "b = ";
        b.print();

        x = y;

    }


    cout << endl;

    return 0;

}

