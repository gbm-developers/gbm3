// header file for matrix template class
// NOTE:  all matrices handled here must be SQUARE 
//        (i.e., # rows = # columns)
// in addition, all DIAGONAL ELEMENTS MUST BE NONZERO
// written by Mike Dinolfo 12/98
// version 1.0

#ifndef __mjdmatrix_h
#define __mjdmatrix_h
#include <iostream>
// generic object (class) definition of matrix:
template <class D> class matrix{
  // NOTE: maxsize determines available memory storage, but
  // actualsize determines the actual size of the stored matrix in use
  // at a particular time.
    int maxsize;  // max number of rows (same as max number of columns)
    int actualsize;  // actual size (rows, or columns) of the stored matrix
    D* data;      // where the data contents of the matrix are stored
    void allocateD()   
    {
        delete[] data;
        data = new D [maxsize*maxsize];
    };
public:
    matrix() 
    {
        maxsize = 5; 
        actualsize = 5;
        data = 0;
        allocateD();
    };                  // private ctor's

    matrix(int newmaxsize) {matrix(newmaxsize,newmaxsize);};
    matrix(int newmaxsize, int newactualsize)  
    { // the only public ctor
        if (newmaxsize <= 0) newmaxsize = 5;
        maxsize = newmaxsize; 
        if ((newactualsize <= newmaxsize)&&(newactualsize>0))
            actualsize = newactualsize;
        else 
            actualsize = newmaxsize;
        // since allocateD() will first call delete[] on data:
        data = 0;
        allocateD();
    };

    ~matrix() 
    { 
        delete[] data; 
    };

    void dumpMatrixValues()
    {
        bool xyz;
        double rv;
        for (int i=0; i < actualsize; i++)
        {
            cout << "i=" << i << ": ";
            for (int j=0; j<actualsize; j++)
            {
                getvalue(i,j,rv,xyz);
                cout << rv << " ";
            } 
            cout << endl;
        }
    };

    void comparetoidentity()  
    {
        int worstdiagonal = 0;
        D maxunitydeviation = 0.0;
        D currentunitydeviation;
        int i = 0;
        for ( i = 0; i < actualsize; i++ )  
        {
            currentunitydeviation = data[i*maxsize+i] - 1.;
            if ( currentunitydeviation < 0.0) currentunitydeviation *= -1.;
            if ( currentunitydeviation > maxunitydeviation )  
            {
                maxunitydeviation = currentunitydeviation;
                worstdiagonal = i;
            }
        }
        int worstoffdiagonalrow = 0;
        int worstoffdiagonalcolumn = 0;
        D maxzerodeviation = 0.0;
        D currentzerodeviation ;
        for ( i = 0; i < actualsize; i++ )  
        {
            for ( int j = 0; j < actualsize; j++ )  
            {
                if ( i == j ) continue;  // we look only at non-diagonal terms
                currentzerodeviation = data[i*maxsize+j];
                if ( currentzerodeviation < 0.0) currentzerodeviation *= -1.0;
                if ( currentzerodeviation > maxzerodeviation )  
                {
                    maxzerodeviation = currentzerodeviation;
                    worstoffdiagonalrow = i;
                    worstoffdiagonalcolumn = j;
                }
            }
        }
        cout << "Worst diagonal value deviation from unity: " 
             << maxunitydeviation << " at row/column " << worstdiagonal << endl;
        cout << "Worst off-diagonal value deviation from zero: " 
             << maxzerodeviation << " at row = " << worstoffdiagonalrow 
             << ", column = " << worstoffdiagonalcolumn << endl;
    }

    void settoproduct(matrix& left, matrix& right)  
    {
        actualsize = left.getactualsize();
        if ( maxsize < left.getactualsize() )   
        {
            maxsize = left.getactualsize();
            allocateD();
        }
        for ( int i = 0; i < actualsize; i++ )
        {
            for ( int j = 0; j < actualsize; j++ )  
            {
                D sum = 0.0;
                D leftvalue, rightvalue;
                bool success;
                for (int c = 0; c < actualsize; c++)  
                {
                    left.getvalue(i,c,leftvalue,success);
                    right.getvalue(c,j,rightvalue,success);
                    sum += leftvalue * rightvalue;
                }
                setvalue(i,j,sum);
            }
        }
    }

    void copymatrix(matrix&  source)  
    {
        actualsize = source.getactualsize();
        if ( maxsize < source.getactualsize() )  
        {
            maxsize = source.getactualsize();
            allocateD();
        }
        for ( int i = 0; i < actualsize; i++ )
        {
            for ( int j = 0; j < actualsize; j++ )  
            {
                D value;
                bool success;
                source.getvalue(i,j,value,success);
                data[i*maxsize+j] = value;
            }
        }
    };

    void setactualsize(int newactualsize) 
    {
        if ( newactualsize > maxsize )
        {
            maxsize = newactualsize ; // * 2;  // wastes memory but saves
                                          // time otherwise required for
                                          // operation new[]
            allocateD();
        }
        if (newactualsize >= 0) actualsize = newactualsize;
    };

    int getactualsize() { return actualsize; };

    void getvalue(int row, int column, D& returnvalue, bool& success)
    {
        if ( (row>=maxsize) || (column>=maxsize) || (row<0) || (column<0) )
        {  
            success = false;
            return;    
        }
        returnvalue = data[ row * maxsize + column ];
        success = true;
    };

    bool setvalue(int row, int column, D newvalue)  
    {
        if ( (row >= maxsize) || (column >= maxsize) || (row<0) || (column<0) )
            return false;
        data[ row * maxsize + column ] = newvalue;
        return true;
    };

    void invert()  
    {
        int i = 0;
        int j = 0;
        int k = 0;
        if (actualsize <= 0) return;  // sanity check
        if (actualsize == 1)
        {
            data[0] = 1.0/data[0];
            return;
        }
        for (i=1; i < actualsize; i++) data[i] /= data[0]; // normalize row 0
        for (i=1; i < actualsize; i++)  
        { 
            for ( j=i; j < actualsize; j++)  
            { // do a column of L
                D sum = 0.0;
                for ( k = 0; k < i; k++)  
                    sum += data[j*maxsize+k] * data[k*maxsize+i];
                data[j*maxsize+i] -= sum;
            }
            if (i == actualsize-1) continue;
            for ( j=i+1; j < actualsize; j++)  
            {  // do a row of U
                D sum = 0.0;
                for ( k = 0; k < i; k++)
                    sum += data[i*maxsize+k]*data[k*maxsize+j];
                data[i*maxsize+j] = (data[i*maxsize+j]-sum) / data[i*maxsize+i];
            }
        }
        for ( i = 0; i < actualsize; i++ )  // invert L
        {
            for ( j = i; j < actualsize; j++ )  
            {
                D x = 1.0;
                if ( i != j ) 
                {
                    x = 0.0;
                    for ( k = i; k < j; k++ ) 
                        x -= data[j*maxsize+k]*data[k*maxsize+i];
                }
                data[j*maxsize+i] = x / data[j*maxsize+j];
            }
        }
        for ( i = 0; i < actualsize; i++ )   // invert U
        {
            for (  j = i; j < actualsize; j++ )  
            {
                if ( i == j ) continue;
                D sum = 0.0;
                for ( k = i; k < j; k++ )
                    sum += data[k*maxsize+j]*( (i==k) ? 1.0 : data[i*maxsize+k] );
                data[i*maxsize+j] = -sum;
            }
        }
        for ( i = 0; i < actualsize; i++ )   // final inversion
        {
            for ( j = 0; j < actualsize; j++ )  
            {
                D sum = 0.0;
                for ( k = ((i>j)?i:j); k < actualsize; k++ )  
                    sum += ((j==k)?1.0:data[j*maxsize+k])*data[k*maxsize+i];
                data[j*maxsize+i] = sum;
            }
        }
    };
};
#endif



