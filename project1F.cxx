#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#define NORMALS

using std::cerr;
using std::endl;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

double norm(double *A) {
    return sqrt((A[0] * A[0]) + (A[1] * A[1]) + (A[2] * A[2]));
}

void normalize(double *A) {
    double _norm = norm(A);
    for (int i = 0; i < 3; i++) {
        A[i] = A[i] / _norm;
    }
}

double dotproduct(double *A, double *B) {
    return ((A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2]));
}

void crossproduct(double *A, double *B, double *AB) {
    AB[0] = (A[1] * B[2]) - (A[2] * B[1]);
    AB[1] = (B[0] * A[2]) - (A[0] * B[2]);
    AB[2] = (A[0] * B[1]) - (A[1] * B[0]);
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Triangle
{
public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double colors[3][3];
    double normals[3][3];
    double shading [3];
};

class Matrix
{
public:
    double          A[4][4];  // A[i][j] means row i, column j
    
    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    
    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
    + ptIn[1]*A[1][0]
    + ptIn[2]*A[2][0]
    + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
    + ptIn[1]*A[1][1]
    + ptIn[2]*A[2][1]
    + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
    + ptIn[1]*A[1][2]
    + ptIn[2]*A[2][2]
    + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
    + ptIn[1]*A[1][3]
    + ptIn[2]*A[2][3]
    + ptIn[3]*A[3][3];
}

class Camera
{
public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    Matrix          finalm;
    
    Matrix DeviceTransform(int sWidth, int sHeight) {
        Matrix m;
        
        m.A[0][0] = sWidth / 2;
        m.A[0][1] = 0;
        m.A[0][2] = 0;
        m.A[0][3] = 0;
        
        m.A[1][0] = 0;
        m.A[1][1] = sHeight / 2;
        m.A[1][2] = 0;
        m.A[1][3] = 0;
        
        m.A[2][0] = 0;
        m.A[2][1] = 0;
        m.A[2][2] = 1;
        m.A[2][3] = 0;
        
        m.A[3][0] = sWidth / 2;
        m.A[3][1] = sHeight / 2;
        m.A[3][2] = 0;
        m.A[3][3] = 1;
        
        return m;
    };
    
    Matrix CameraTransform(double *u, double *v, double *w, double *O) {
        Matrix m;
        
        m.A[0][0] = u[0];
        m.A[0][1] = v[0];
        m.A[0][2] = w[0];
        m.A[0][3] = 0;
        
        m.A[1][0] = u[1];
        m.A[1][1] = v[1];
        m.A[1][2] = w[1];
        m.A[1][3] = 0;
        
        m.A[2][0] = u[2];
        m.A[2][1] = v[2];
        m.A[2][2] = w[2];
        m.A[2][3] = 0;
        
        double t[3];
        for (int i = 0; i < 3; i++) {
            t[i] = 0 - O[i];
        }
        
        m.A[3][0] = dotproduct(u, t);
        m.A[3][1] = dotproduct(v, t);
        m.A[3][2] = dotproduct(w, t);
        m.A[3][3] = 1;
        
        return m;
    };
    
    Matrix ViewTransform(void) {
        Matrix m;
        
        m.A[0][0] = 1 / tan(angle/2);
        m.A[0][1] = 0;
        m.A[0][2] = 0;
        m.A[0][3] = 0;
        
        m.A[1][0] = 0;
        m.A[1][1] = 1 / tan(angle/2);
        m.A[1][2] = 0;
        m.A[1][3] = 0;
        
        m.A[2][0] = 0;
        m.A[2][1] = 0;
        m.A[2][2] = (far + near) / (far - near);
        m.A[2][3] = -1;
        
        m.A[3][0] = 0;
        m.A[3][1] = 0;
        m.A[3][2] = 2 * (far * near) / (far - near);
        m.A[3][3] = 0;
        
        return m;
    };
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    
    double u[3]; double v[3]; double w[3];
    for (int i = 0; i < 3; i++) {
        w[i] = c.position[i] - c.focus[i];
    }
    normalize(w);
    crossproduct(c.up, w, u);
    normalize(u);
    crossproduct(w, u, v);
    normalize(v);
    Matrix dt = c.DeviceTransform(1000, 1000);
    Matrix ct = c.CameraTransform(u, v, w, c.position);
    Matrix vt = c.ViewTransform();
    c.finalm = c.finalm.ComposeMatrices(ct, vt);
    c.finalm = c.finalm.ComposeMatrices(c.finalm, dt);
    //c.finalm.Print(cout);
    //cout << endl;
    return c;
}

class Screen
{
public:
    unsigned char   *buffer;
    int width, height;
    double *Zbuffer;
    Camera camera;
};

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif
        
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
            { 0, 0, 91 },
            { 0, 255, 255 },
            { 0, 128, 0 },
            { 255, 255, 0 },
            { 255, 96, 0 },
            { 107, 0, 0 },
            { 224, 76, 76 }
        };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }
    
    return tris;
}

double CalculateShading(double *normal, Triangle *t, Camera cam, int idx) {
    LightingParameters LP;
    double diffuse = fabs(dotproduct(LP.lightDir, normal));
    double R[3];
    R[0] = (dotproduct(LP.lightDir, normal) * 2.0) * normal[0] - LP.lightDir[0];
    R[1] = (dotproduct(LP.lightDir, normal) * 2.0) * normal[1] - LP.lightDir[1];
    R[2] = (dotproduct(LP.lightDir, normal) * 2.0) * normal[2] - LP.lightDir[2];
    double viewDirection[3];
    viewDirection[0] = cam.position[0] - t->X[idx];
    viewDirection[1] = cam.position[1] - t->Y[idx];
    viewDirection[2] = cam.position[2] - t->Z[idx];
    normalize(viewDirection);
    double specular = fmax(0, pow((dotproduct(R, viewDirection)), LP.alpha));
    double shadingAmount = LP.Ka + LP.Kd * diffuse + LP.Ks * specular;
    return shadingAmount;
}

double LERP(double tt, double b, double e) {
    return tt * (e - b) + b;
}

void SetPixels(Triangle *t, Screen *s, int scanline, double leftEnd, double rightEnd, double tt) {
    double r1 = LERP(tt, t->colors[1][0], t->colors[0][0]);
    double r2 = LERP(tt, t->colors[2][0], t->colors[0][0]);
    double g1 = LERP(tt, t->colors[1][1], t->colors[0][1]);
    double g2 = LERP(tt, t->colors[2][1] ,t->colors[0][1]);
    double b1 = LERP(tt, t->colors[1][2] ,t->colors[0][2]);
    double b2 = LERP(tt, t->colors[2][2] ,t->colors[0][2]);
    double z1 = LERP(tt, t->Z[1], t->Z[0]);
    double z2 = LERP(tt, t->Z[2], t->Z[0]);
    double s1 = LERP(tt, t->shading[1], t->shading[0]);
    double s2 = LERP(tt, t->shading[2], t->shading[0]);
    for (int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++) {
        if (c == s->width) {return;}
        int index = 3 * (c + scanline * s->width);
        if (index < 0 or index > (3 * s->width * s->height)) { break; }
        tt = (c - leftEnd) / (rightEnd - leftEnd);
        if (s->Zbuffer[index] <= LERP(tt, z1, z2)) {
            s->Zbuffer[index] = LERP(tt, z1, z2);
            double shading = LERP(tt, s1, s2);
            s->buffer[index] = fmin(255.0, ceil_441(shading * LERP(tt, r1, r2) * 255.0));
            s->buffer[index+1] = fmin(255.0, ceil_441(shading * LERP(tt, g1, g2) * 255.0));
            s->buffer[index+2] = fmin(255.0, ceil_441(shading * LERP(tt, b1, b2) * 255.0));
        }
    }
}

void TriangleSplit(Triangle *t, Triangle *top, Triangle *bot) {
    int iTop = -1; int iBot = -1; int iMid = -1; int iLeft = -1; int iRight = -1;
    int size = 3* sizeof(double);
    double x0 = t->X[0]; double x1 = t->X[1]; double x2 = t->X[2];
    double y0 = t->Y[0]; double y1 = t->Y[1]; double y2 = t->Y[2];
    
    if (y0 == y1 or y1 == y2 or y0 == y2) {
        if (y0 == y1) {
            if (y2 > y0) {
                iTop = 2;
            }
            else {
                iBot = 2;
            }
            if (x0 > x1) {
                iLeft = 1; iRight = 0;
            }
            else {
                iLeft = 0; iRight = 1;
            }
        }
        else if (y0 == y2) {
            if (y1 > y0) {
                iTop = 1;
            }
            else {
                iBot = 1;
            }
            if (x0 > x2) {
                iLeft = 2; iRight = 0;
            }
            else {
                iLeft = 0; iRight = 2;
            }
        }
        else if (y1 == y2) {
            if (y0 > y1) {
                iTop = 0;
            }
            else {
                iBot = 0;
            }
            if (x1 > x2) {
                iLeft = 2; iRight = 1;
            }
            else {
                iLeft = 1; iRight = 2;
            }
        }
        if (iTop == -1) {
            top->X[0] = -1;
            bot->X[0] = t->X[iBot]; bot->X[1] = t->X[iLeft]; bot->X[2] = t->X[iRight];
            bot->Y[0] = t->Y[iBot]; bot->Y[1] = t->Y[iLeft]; bot->Y[2] = t->Y[iRight];
            bot->Z[0] = t->Z[iBot]; bot->Z[1] = t->Z[iLeft]; bot->Z[2] = t->Z[iRight];
            bot->shading[0] = t->shading[iBot]; bot->shading[1] = t->shading[iLeft]; bot->shading[2] = t->shading[iRight];
            memcpy(&bot->colors[0], &t->colors[iBot], size);
            memcpy(&bot->colors[1], &t->colors[iLeft], size);
            memcpy(&bot->colors[2], &t->colors[iRight], size);
        }
        else {
            bot->X[0] = -1;
            top->X[0] = t->X[iTop]; top->X[1] = t->X[iLeft]; top->X[2] = t->X[iRight];
            top->Y[0] = t->Y[iTop]; top->Y[1] = t->Y[iLeft]; top->Y[2] = t->Y[iRight];
            top->Z[0] = t->Z[iTop]; top->Z[1] = t->Z[iLeft]; top->Z[2] = t->Z[iRight];
            top->shading[0] = t->shading[iBot]; top->shading[1] = t->shading[iLeft]; top->shading[2] = t->shading[iRight];
            memcpy(&top->colors[0], &t->colors[iTop], size);
            memcpy(&top->colors[1], &t->colors[iLeft], size);
            memcpy(&top->colors[2], &t->colors[iRight], size);
        }
        return;
    }
    
    if ( (y0 < y1 and y0 > y2) or (y0 < y2 and y0 > y1) ) {
        iMid = 0;
        if (y1 > y2) {
            iTop = 1; iBot = 2;
        }
        else {
            iTop = 2; iBot = 1;
        }
    }
    else if ( (y1 < y2 and y1 > y0) or (y1 < y0 and y1 > y2) ) {
        iMid = 1;
        if (y0 > y2) {
            iTop = 0; iBot = 2;
        }
        else {
            iTop = 2; iBot = 0;
        }
    }
    else if ( (y2 < y0 and y2 > y1) or (y2 < y1 and y2 > y0) ) {
        iMid = 2;
        if (y0 > y1) {
            iTop = 0; iBot = 1;
        }
        else {
            iTop = 1; iBot = 0;
        }
    }
    
    if (t->X[iTop] < t->X[iBot]) {
        if (t->X[iTop] > t->X[iMid]) {
            iLeft = iMid;
        }
        else {
            iLeft = iTop;
        }
        if (t->X[iMid] < t->X[iBot]) {
            iRight = iBot;
        }
        else {
            iRight = iMid;
        }
    }
    else {
        if (t->X[iBot] < t->X[iMid]) {
            iLeft = iBot;
        }
        else {
            iLeft = iMid;
        }
        if (t->X[iMid] < t->X[iTop]) {
            iRight = iTop;
        }
        else {
            iRight = iMid;
        }
    }
    
    double slope = (t->Y[iTop] - t->Y[iBot]) / (t->X[iTop] - t->X[iBot]);
    double b = t->Y[iTop] - (slope * t->X[iTop]);
    double rem = (t->Y[iMid] - b ) / slope;
    double tt = (t->Y[iMid] - t->Y[iBot]) / (t->Y[iTop] - t->Y[iBot]);
    
    if (isnan(slope) == true or isinf(slope) == true) {
        rem = t->X[iBot];
    }
    
    if(t->X[iMid] >= rem) {
        top->X[1] = rem; top->X[2] = t->X[iMid];
        top->Y[1] = t->Y[iMid]; top->Y[2] = t->Y[iMid];
        top->Z[1] = LERP(tt, t->Z[iBot], t->Z[iTop]); top->Z[2] = t->Z[iMid];
        bot->X[1] = rem; bot->X[2] = t->X[iMid];
        bot->Y[1] = t->Y[iMid]; bot->Y[2] = t->Y[iMid];
        bot->Z[1] = LERP(tt, t->Z[iBot], t->Z[iTop]); bot->Z[2] = t->Z[iMid];
        top->colors[1][0] = LERP(tt, t->colors[iBot][0], t->colors[iTop][0]);
        top->colors[1][1] = LERP(tt, t->colors[iBot][1], t->colors[iTop][1]);
        top->colors[1][2] = LERP(tt, t->colors[iBot][2], t->colors[iTop][2]);
        bot->colors[1][0] = LERP(tt, t->colors[iBot][0], t->colors[iTop][0]);
        bot->colors[1][1] = LERP(tt, t->colors[iBot][1], t->colors[iTop][1]);
        bot->colors[1][2] = LERP(tt, t->colors[iBot][2], t->colors[iTop][2]);
        bot->shading[2] = t->shading[iMid];
        bot->shading[1] = LERP(tt, t->shading[iBot], t->shading[iTop]);
        top->shading[2] = t->shading[iMid];
        top->shading[1] = LERP(tt, t->shading[iBot], t->shading[iTop]);
        memcpy(&top->colors[2], &t->colors[iMid], size);
        memcpy(&bot->colors[2], &t->colors[iMid], size);
    }
    else {
        top->X[1] = t->X[iMid]; top->X[2] = rem;
        top->Y[1] = t->Y[iMid]; top->Y[2] = t->Y[iMid];
        top->Z[1] = t->Z[iMid]; top->Z[2] = LERP(tt, t->Z[iBot], t->Z[iTop]);
        bot->X[1] = t->X[iMid]; bot->X[2] = rem;
        bot->Y[1] = t->Y[iMid]; bot->Y[2] = t->Y[iMid];
        bot->Z[1] = t->Z[iMid]; bot->Z[2] = LERP(tt, t->Z[iBot], t->Z[iTop]);
        top->colors[2][0] = LERP(tt, t->colors[iBot][0], t->colors[iTop][0]);
        top->colors[2][1] = LERP(tt, t->colors[iBot][1], t->colors[iTop][1]);
        top->colors[2][2] = LERP(tt, t->colors[iBot][2], t->colors[iTop][2]);
        bot->colors[2][0] = LERP(tt, t->colors[iBot][0], t->colors[iTop][0]);
        bot->colors[2][1] = LERP(tt, t->colors[iBot][1], t->colors[iTop][1]);
        bot->colors[2][2] = LERP(tt, t->colors[iBot][2], t->colors[iTop][2]);
        bot->shading[1] = t->shading[iMid];
        bot->shading[2] = LERP(tt, t->shading[iBot], t->shading[iTop]);
        top->shading[1] = t->shading[iMid];
        top->shading[2] = LERP(tt, t->shading[iBot], t->shading[iTop]);
        memcpy(&top->colors[1], &t->colors[iMid], size);
        memcpy(&bot->colors[1], &t->colors[iMid], size);
    }
    bot->X[0] = t->X[iBot];
    bot->Y[0] = t->Y[iBot];
    bot->Z[0] = t->Z[iBot];
    top->X[0] = t->X[iTop];
    top->Y[0] = t->Y[iTop];
    top->Z[0] = t->Z[iTop];
    bot->shading[0] = t->shading[iBot];
    top->shading[0] = t->shading[iTop];
    memcpy(&bot->colors[0], &t->colors[iBot], size);
    memcpy(&top->colors[0], &t->colors[iTop], size);
    return;
}

void RasterizeGoingDownTriangle(Triangle *t, Screen *s) {
    double x0 = t->X[0]; double x1 = t->X[1]; double x2 = t->X[2];
    double y0 = t->Y[0]; double y1 = t->Y[1]; double y2 = t->Y[2];
    double slope1 = (x1 - x0) / (y1 - y0);
    double slope2 = (x2 - x0) / (y2 - y0);
    double rowMin = ceil_441(y0);
    double rowMax = floor_441(y1);
    double leftEnd = x0 + slope1 * (rowMin - y0);
    double rightEnd = x0 + slope2 * (rowMin - y0);
    
    for (int r = rowMin; r <= rowMax; r++) {
        double tt = (y1 - r) / (y1 - y0);
        SetPixels(t, s, r, leftEnd, rightEnd, tt);
        leftEnd += slope1;
        rightEnd += slope2;
    }
}

void RasterizeGoingUpTriangle(Triangle *t, Screen *s) {
    double x0 = t->X[0]; double x1 = t->X[1]; double x2 = t->X[2];
    double y0 = t->Y[0]; double y1 = t->Y[1]; double y2 = t->Y[2];
    double slope1 = (x0 - x1) / (y0 - y1);
    double slope2 = (x0 - x2) / (y0 - y2);
    double rowMin = ceil_441(y1);
    double rowMax = floor_441(y0);
    double leftEnd = x1 + slope1 * (rowMin - y1);
    double rightEnd = x2 + slope2 * (rowMin - y1);
    
    for (int r = rowMin; r <= rowMax; r++) {
        double tt = (r - y1) / (y0 - y1);
        SetPixels(t, s, r, leftEnd, rightEnd, tt);
        leftEnd += slope1;
        rightEnd += slope2;
    }
}

void ApplyMatrix(Triangle *t, Camera cam) {
    double ov1[4]; double ov2[4]; double ov3[4];
    double nv1[4]; double nv2[4]; double nv3[4];
    
    for (int i = 0; i < 3; i++) {
        if (i == 0) {
            ov1[0] = t->X[i];
            ov1[1] = t->Y[i];
            ov1[2] = t->Z[i];
            ov1[3] = 1;
        }
        else if (i == 1) {
            ov2[0] = t->X[i];
            ov2[1] = t->Y[i];
            ov2[2] = t->Z[i];
            ov2[3] = 1;
        }
        else {
            ov3[0] = t->X[i];
            ov3[1] = t->Y[i];
            ov3[2] = t->Z[i];
            ov3[3] = 1;
        }
    }
    
    cam.finalm.TransformPoint(ov1, nv1);
    cam.finalm.TransformPoint(ov2, nv2);
    cam.finalm.TransformPoint(ov3, nv3);
    
    t->X[0] = nv1[0] / nv1[3]; t->X[1] = nv2[0] / nv2[3]; t->X[2] = nv3[0] / nv3[3];
    t->Y[0] = nv1[1] / nv1[3]; t->Y[1] = nv2[1] / nv2[3]; t->Y[2] = nv3[1] / nv3[3];
    t->Z[0] = nv1[2] / nv1[3]; t->Z[1] = nv2[2] / nv2[3]; t->Z[2] = nv3[2] / nv3[3];
}

void RasterizeArbitraryTriangle(Triangle *t, Screen *s, Camera cam) {
    t->shading[0] = CalculateShading(t->normals[0], t, cam, 0);
    t->shading[1] = CalculateShading(t->normals[1], t, cam, 1);
    t->shading[2] = CalculateShading(t->normals[2], t, cam, 2);
    ApplyMatrix(t, cam);
    Triangle top, bot;
    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {
            top.colors[c1][c2] = t->colors[c1][c2];
            bot.colors[c1][c2] = t->colors[c1][c2];
        }
    }
    TriangleSplit(t, &top, &bot);
    if (top.X[0] != -1) {
        RasterizeGoingUpTriangle(&top, s);
        }
    if (bot.X[0] != -1) {
        RasterizeGoingDownTriangle(&bot, s);
        }
}

void MakeImage(char const* name, std::vector<Triangle> triangles, int frame, int nframes) {
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer =
    (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1000*1000;
    for (int i = 0 ; i < npixels*3 ; i++)
        buffer[i] = 0;
    
    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;
    
    double *Zbuffer = new double [3 * screen.width * screen.height];
    for (int j = 0; j < npixels*3; j++)
        Zbuffer[j] = -1;
    screen.Zbuffer = Zbuffer;
    
    Camera cam = GetCamera(frame, nframes);
    screen.camera = cam;
    
    for (std::vector<Triangle>::iterator it = triangles.begin(); it != triangles.end(); it++) {
        RasterizeArbitraryTriangle(&(*it), &screen, cam);
    }
    WriteImage(image, name);
}

int main() {
    std::vector<Triangle> triangles = GetTriangles();
    MakeImage("frame000", triangles, 0, 1000);
    //MakeImage("frame250", triangles, 250, 1000);
    //MakeImage("frame500", triangles, 500, 1000);
    //MakeImage("frame750", triangles, 750, 1000);
    return 1;
}
