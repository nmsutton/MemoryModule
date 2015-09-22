: $Id: matrix.mod,v 1.32 2003/09/05 18:52:54 billl Exp $

:* COMMENT
COMMENT
NB: only minimal error checking done
NB: no dynamic allocation - eg m1.transpose(m1) will give a wrong result 
NB: matrix and vec sizes must be correct before using: use .resize()

================            USAGE               ================
objref mat
mat = new Vector(rows*cols)
mat.mprintf(M,N)     // print out as M rows and N columns
mat2.transpose(mat)  // transpose of matrix
y.mmult(mat,x)       // y = mat*x
y.spmult(pre,post,mat,x) // y = mat*x using "sparse matrix"
w.spget(pre,post,row,col) // ie pre,post,post,pre!!
spidget(pre,post,prid,poid ...) // match vals in 4 vecs to 4 vals
wt.mkspcp/chkspcp(pre,post) // copy the indices into integer arrays
mat.outprod(x,y)     // mat = outer product of vectors x and y
mat.mget(i,j,cols)   // i=row#; j=col#
mat.mset(i,j,cols,val)
y.mrow(mat,i,cols)
y.mcol(mat,j,cols)
y.msetrow(mat,i,cols)
y.msetcol(mat,j,cols)
================================================================
ENDCOMMENT

NEURON {
  SUFFIX nothing
}
 
PARAMETER {
  MATRIX_INSTALLED=0
}

:* mat.outprod(x,y) // mat = outer product of vectors x and y
VERBATIM
static double outprod(void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z;
  /* this will be the outer product */
  nx = vector_instance_px(vv, &x);
	
  /* these are the two vectors that make it up */
  ny = vector_arg_px(1, &y); // will be number of columns
  nz = vector_arg_px(2, &z); // will be number of rows
  if (nx != ny*nz) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<ny;i++) {
    for (j=0;j<nz;j++) {
      x[i*nz+j] = y[i]*z[j];
    }
  }
  return nx;
}
ENDVERBATIM
 
:* mmult
VERBATIM
static double mmult(void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z;
  /* x will be the product of matrix y and vec z */
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  nz = vector_arg_px(2, &z);
  if (ny != nx*nz) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<nx;i++) {
    x[i] = 0.;
    for (j=0;j<nz;j++) {
      x[i] += y[i*nz+j]*z[j];
    }
  }
  return nx;
}
ENDVERBATIM
 
:* ST[PO].spltp(pr,po,wt,ST[PRE])
VERBATIM
static double spltp(void* vv) {
  int ii, jj, nstpr, nstpo, nw, npr, npo, flag, cnt;
  double *stpr, *stpo, *w, *pr, *po;
  extern double hoc_call_func(Symbol*, int narg);

  char func[4] = "ltp";
  Symbol* s = hoc_lookup(func);
  if (! s) { hoc_execerror("Can't find ltp() func", 0); }
  nstpo = vector_instance_px(vv, &stpo);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nw = vector_arg_px(3, &w);
  nstpr = vector_arg_px(4, &stpr);
  for (ii=0,jj=0,cnt=0;ii<nstpo;ii++) {
    if (stpo[ii]==1.0) { /* connections to these will be changed */ 
      for (;po[jj]<ii;jj++) ; /* move forward till find a po */
      for (;po[jj]==ii;jj++) { /* move through these po's */
	if (stpr[(int)pr[jj]]==1.) { /*  did the presyn spike? */
	  cnt++; hoc_pushx(1.0);
	} else { 
	  cnt--; hoc_pushx(-1.0);
	}
        hoc_pushx(w[jj]);
        w[jj]=hoc_call_func(s, 2);
      }
    }
  }
  return cnt;
}
ENDVERBATIM
 
VERBATIM
/* Maintain a parallel vector of ints to avoid the slowness of repeated casts in spmult */
static int *pr_int;
static int *po_int;
static int cpfl=0;
ENDVERBATIM

:* wt.mkspcp(pr,po)
VERBATIM
static double mkspcp(void* vv) {
  int j, nw, npr, npo;
  double *w, *pr, *po;
  if (! ifarg(1)) { 
    cpfl=0; 
    if (po_int!=NULL) free(po_int); 
    if (pr_int!=NULL) free(pr_int);
    po_int=(int *)NULL; pr_int=(int *)NULL; 
    return 0;
  }
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  pr_int=(int *)ecalloc(nw, sizeof(int));
  po_int=(int *)ecalloc(nw, sizeof(int));
  for (j=0;j<nw;j++) {
    po_int[j]=(int)po[j];
    pr_int[j]=(int)pr[j];
  }
  cpfl=nw;
  return cpfl;
}
ENDVERBATIM

:* wt.chkspcp(pr,po)
VERBATIM
static double chkspcp(void* vv) {
  int j, nw, npr, npo, flag;
  double *w, *pr, *po;
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  flag=1;
  if (po_int==NULL || pr_int==NULL) { cpfl=0; return 0; }
  if (cpfl!=nw) { flag=0;
  } else for (j=0;j<nw;j++) {
    if (po_int[j]!=(int)po[j] || pr_int[j]!=(int)pr[j]) {flag=0; continue;}
  }
  if (flag==0) {
    cpfl=0; free(po_int); free(pr_int); 
    po_int=(int *)NULL; pr_int=(int *)NULL; 
  }
  return flag;
}
ENDVERBATIM

:* y.spmult(pr,po,wt,x[,flag])
: y=W*x, y will be the product of matrix w with pre/post indices and vec x
: optional flag (5th arg present) - do not clear dest vector initially
VERBATIM
static double spmult(void* vv) {
  int i, j, nx, ny, nw, npr, npo, flag;
  double *x, *y, *w, *pr, *po, xx;
  ny = vector_instance_px(vv, &y);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nw = vector_arg_px(3, &w);
  nx = vector_arg_px(4, &x);
  if (ifarg(5)) {flag=1;} else {flag=0;}
  if (nw!=npr || nw!=npo) {
    hoc_execerror("Sparse mat must have 3 identical size vecs for pre/post/wt", 0); 
  }
  if (flag==0) for (i=0;i<ny;i++) y[i] = 0.; // clear dest vec
  if (cpfl==0) {
    for (j=0;j<nw;j++) y[(int)po[j]] += (x[(int)pr[j]]*w[j]);
  } else if (cpfl!=nw) { hoc_execerror("cpfl!=nw in spmult", 0); } else {
    for (j=0;j<nw;j++) if (x[pr_int[j]]!=0) { y[po_int[j]] += ((x[pr_int[j]])*w[j]); }
  }
  return nx;
}
ENDVERBATIM
 

:* wt.spget(pr,po,row,col) returns weight value
VERBATIM
static double spget(void* vv) {
  int j, nw, npr, npo;
  double *w, *pr, *po, row, col;
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  row = *getarg(3);
  col = *getarg(4);
  for (j=0;j<nw;j++) if (row==po[j]&&col==pr[j]) break;
  if (j==nw) return 0.; else return w[j];
}
ENDVERBATIM

:* spidget(prv,pov,pridv,poidv,pr,po,prid,poid) returns an index
FUNCTION spidget() { : gaussian distribution around 0
VERBATIM
{
  int j, npr, npo, nprid, npoid;
  double *pr, *po, *prid, *poid, pri, poi, pridi, poidi;
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nprid = vector_arg_px(3, &prid);
  npoid = vector_arg_px(4, &poid);
  pri= *getarg(5);
  poi= *getarg(6);
  pridi= *getarg(7);
  poidi= *getarg(8);
  for (j=0;j<npr;j++) { 
    if (poi==po[j]&&pri==pr[j]&&pridi==prid[j]&&poidi==poid[j]) break;
  }
  if (j==npr) _lspidget=-1.0; else _lspidget=(double)j;
}
ENDVERBATIM
}

:* transpose
VERBATIM
static double transpose(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  /* x will be the transpose of matrix y */
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  rows = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (ny != nx) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<rows;i++) {
    for (j=0;j<cols;j++) {
      x[j*rows+i] = y[i*cols+j];
    }
  }
  return nx;
}
ENDVERBATIM
 
:* mprintf
VERBATIM
static double mprintf(void* vv) {
  int i, j, nx, rows, cols;
  double *x;
  /* x will be printed out */
  nx = vector_instance_px(vv, &x);
  rows = (int)*getarg(1);
  cols = (int)*getarg(2);
  if (nx != rows*cols) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<rows;i++) {
    for (j=0;j<cols;j++) {
      printf("%g\t",x[i*cols+j]);
    }
    printf("\n");
  }
  return nx;
}
ENDVERBATIM
 
:* mget(i,j,cols)
VERBATIM
static double mget(void* vv) {
  int i, j, nx, rows, cols;
  double *x;
  nx = vector_instance_px(vv, &x);
  i = (int)*getarg(1);
  j = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (i*cols+j >= nx) {
    hoc_execerror("Indices out of bounds", 0);
  }
  return x[i*cols+j];
}
ENDVERBATIM
 
:* mrow(mat,i,cols)
VERBATIM
static double mrow(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  i = (int)*getarg(2);
  cols = (int)*getarg(3); rows=ny/cols;
  if (cols!=nx) { 
    nx=vector_buffer_size(vv);
    if (cols<=nx) {
      vector_resize(vv, cols); nx=cols; 
    } else {
      printf("%d > %d :: \n",cols,nx);
      hoc_execerror("Vector max capacity too small in mrow", 0);
    }
  }
  if (i>=rows) { hoc_execerror("Indices out of bounds", 0); }
  for (j=0;j<nx;j++) { x[j] = y[i*cols+j]; }
  return nx;
}
ENDVERBATIM
 
:* mcol(mat,j,cols)
VERBATIM
static double mcol(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  j = (int)*getarg(2);
  cols = (int)*getarg(3); rows=ny/cols;
  if (rows!=nx) { 
    nx=vector_buffer_size(vv);
    if (rows<=nx) {
      vector_resize(vv, rows); nx=rows; 
    } else {
      printf("%d > %d :: ",rows,nx);
      hoc_execerror("Vector max capacity too small in mcol ", 0);
    }
  }
  if (j>=cols) { hoc_execerror("Indices out of bounds", 0); }
  for (i=0;i<nx;i++) { x[i] = y[i*cols+j]; }
  return nx;
}
ENDVERBATIM
 
:* msetrow(mat,i,cols)
VERBATIM
static double msetrow(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  i = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (cols!=nx || i>=ny/cols) {
    hoc_execerror("Indices out of bounds", 0);
  }
  for (j=0;j<nx;j++) { y[i*cols+j] = x[j]; }
  return nx;
}
ENDVERBATIM

:* msetcol(mat,j,cols)
VERBATIM
static double msetcol(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  j = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (cols!=ny/nx || j>=cols) {
    hoc_execerror("Indices out of bounds", 0);
  }
  for (i=0;i<nx;i++) { y[i*cols+j] = x[i]; }
  return nx;
}
ENDVERBATIM
 
:* mset(i,j,cols,val)
VERBATIM
static double mset(void* vv) {
  int i, j, nx, rows, cols;
  double *x, val;
  nx = vector_instance_px(vv, &x);
  i = (int)*getarg(1);
  j = (int)*getarg(2);
  cols = (int)*getarg(3);
  val = *getarg(4);
  if (i*cols+j >= nx) {
    hoc_execerror("Indices out of bounds", 0);
  }
  return (x[i*cols+j]=val);
}
ENDVERBATIM
 
:* PROCEDURE install_matrix()
PROCEDURE install_matrix() {
  MATRIX_INSTALLED=1
VERBATIM
  /* the list of additional methods */
  install_vector_method("outprod", outprod);
  install_vector_method("mmult", mmult);
  install_vector_method("spmult", spmult);
  install_vector_method("spget", spget);
  install_vector_method("mkspcp", mkspcp);
  install_vector_method("chkspcp", chkspcp);
  install_vector_method("spltp", spltp);
  install_vector_method("transpose", transpose);
  install_vector_method("mprintf", mprintf);
  install_vector_method("mget", mget);
  install_vector_method("mset", mset);
  install_vector_method("mrow", mrow);
  install_vector_method("mcol", mcol);
  install_vector_method("msetrow", msetrow);
  install_vector_method("msetcol", msetcol);
ENDVERBATIM
}
