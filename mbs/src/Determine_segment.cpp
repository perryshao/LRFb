// MATLAB MEX gateway for projected maximal blurred segment recognition.
//

#include <mex.h>
#include <matrix.h>
#include "unit_reco.cpp"

using namespace std;

CMTMemoryPool<CMemoryPool<Reco>, CCriticalSection> *Reco::s_pool = NULL;

#define Fori(x) for (int i = 0; i < (x); i++)
#define Forj(x) for (int j = 0; j < (x); j++)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *P, *width;
    double *bSegment, *eSegment;
    int m, pos;

    /* parse input arguments*/
    P = mxGetPr(prhs[0]);     //P
    m = mxGetM(prhs[0]);      //rows for P
    width = mxGetPr(prhs[1]); //width
    point P_i;

    /* create output arguments*/
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //bSegment
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); //eSegment
    bSegment = mxGetPr(plhs[0]);
    eSegment = mxGetPr(plhs[1]);

    /* define the parameters*/
    Reco::NewPool();
    Reco *M[4];
    M[0] = new RecoLR();
    M[1] = new RecoRL();
    M[2] = new RecoTB();
    M[3] = new RecoBT();
    Q curEp;
    Q max_thickness((Z)width[0], 1);
    pos = 0;
    /* initialize  */
    int typemax = 0;
    int indmax = 0;

    Fori(4)
    {
        M[i]->Init();
        curEp = 0;
        int k = 0;
        while ((curEp < max_thickness) && (pos + k < m))
        {
            P_i.x = (Z)P[0 * m + (pos + k)];
            P_i.y = (Z)P[1 * m + (pos + k)];
            M[i]->Insert(curEp, P_i);
            if (curEp < max_thickness)
            {
                if ((indmax == 0) && (M[i]->valid()))
                {
                    indmax = 1;
                    typemax = i;
                }
                if ((indmax > 0) && (k >= indmax))
                {
                    indmax = k + 1;
                    typemax = i;
                }
                k++;
            }
        }
    }

    /*M.Init();
    Fori(indmax)
    {
		P_i.x=P[0*m+(pos+i)]; P_i.y=P[1*m+(pos+i)];
		M.Insert(curEp, P_i);
	}*/
    Fori(4) delete M[i];
    Reco::DeletePool();

    bSegment[0] = pos;
    bSegment[0]++; // keep consistent with Matlab
    pos = pos + indmax - 1;
    eSegment[0] = pos;
    eSegment[0]++; // keep consistent with Matlab
}
