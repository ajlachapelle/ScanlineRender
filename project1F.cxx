#include <iostream>
#include <math.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

using std::cerr;
using std::endl;

#define NORMALS

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}


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

double*
scaleVector(double s, double* v, double* vScaled)
{
  for (int i = 0; i < 3; ++i)
    vScaled[i] = s * v[i];
  return vScaled;
}

double*
vectorDifference(double* v1, double* v2, double* difference)
{
  for (int i = 0; i < 3; ++i)
    difference[i] = v1[i] - v2[i];
  return difference;
}

double
dotProduct(double* v1, double* v2)
{
  double dotProduct = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; // Could be implemented as a for loop over a += operation to account for any vector dimensions
  return dotProduct;
}

double*
crossProduct(double* v1, double* v2, double* product)
{
  product[0] = v1[1]*v2[2] - v1[2]*v2[1];
  product[1] = v1[2]*v2[0] - v1[0]*v2[2];
  product[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return product;
}

double*
normalize(double* v, double* vUnit)
{
  double norm = sqrt(dotProduct(v,v));
  vUnit[0] = v[0]/norm;
  vUnit[1] = v[1]/norm;
  vUnit[2] = v[2]/norm;
  return vUnit;
}

// Returns the interpolated field value of a point between two given endpoints a and b
double interpolate(double a, double Fa, double b, double Fb, double x)
{
  double t = (x - a)/(b - a);
  double Fx = Fa + t*(Fb - Fa);
  return Fx;
}


class Triangle
{
  public:
      double          X[3];
      double          Y[3];
      double          Z[3];
      double  colors[3][3];
      double normals[3][3];
      double    shading[3];

      void ySortVertices(void);

  private:
    void swapVertices(int n, int m);

  // would some methods for transforming the triangle in place be helpful?
};

// Bubble Sort vertices according to Y-coordinate
void
Triangle::ySortVertices(void)
{
  bool sorted = false; // Flag for tracking whether a single iteration of the sort has swapped any elements
  for (int m = 2; m > 0 && !sorted; --m) // If no swaps are made during a single iteration, the matrix is already sorted
  {
    sorted = true;
    for (int n = 0; n < m; ++n)
    {
      if (Y[n] > Y[n+1])
      {
        swapVertices(n,n+1);
        sorted = false;
      }
    }
  }
}

// Helper method for swapping the X, Y, Z, and color values of two vertices
void
Triangle::swapVertices(int n, int m)
{
// With the extra data for each vertex that must be swapped, is it more efficient to just store vertex order seperately?
  double tempX = X[n];
  double tempY = Y[n];
  double tempZ = Z[n];
  double tempColor[3] = {colors[n][0], colors[n][1], colors[n][2]};
  //for (int i = 0; i <= 2; ++i)
    //tempColor[i] = colors[n][i];
  double tempNormal[3] = {normals[n][0], normals[n][1], normals[n][2]};
  double tempShading = shading[n];

  X[n] = X[m];
  Y[n] = Y[m];
  Z[n] = Z[m];
  shading[n] = shading[m];
  for (int i = 0; i <= 2; ++i)
  {
    colors[n][i] = colors[m][i];
    normals[n][i] = normals[m][i];
  }

  X[m] = tempX;
  Y[m] = tempY;
  Z[m] = tempZ;
  shading[m] = tempShading;
  for (int i = 0; i <= 2; ++i)
  {
    colors[m][i] = tempColor[i];
    normals[m][i] = tempNormal[i];
  }
}

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    Matrix(void);

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

Matrix::Matrix(void)
{
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      A[i][j] = 0;
}

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

    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(int width, int height);
};

Matrix
Camera::ViewTransform(void)
{
  Matrix transform;

  transform.A[0][0] = transform.A[1][1] = 1/tan(angle/2);
  transform.A[2][2] = (far+near)/(far-near);
  transform.A[2][3] = -1;
  transform.A[3][2] = 2*far*near/(far-near);

  return transform;
}

Matrix
Camera::CameraTransform(void)
{
  Matrix transform;

  double adjustedFocus[3];
  vectorDifference(position, focus, adjustedFocus);
  double u[3];
  crossProduct(up, adjustedFocus, u);
  normalize(u, u);
  double v[3];
  crossProduct(adjustedFocus, u, v);
  normalize(v, v);
  double w[3] = {adjustedFocus[0], adjustedFocus[1], adjustedFocus[2]};
  normalize(w, w);

  double t[3] = {0};
  vectorDifference(t, position, t);

  for(int i = 0; i <= 2; ++i)
  {
    transform.A[i][0] = u[i];
    transform.A[i][1] = v[i];
    transform.A[i][2] = w[i];
    //transform.A[i][3] = 0; (implicitly)
  }
  transform.A[3][0] = dotProduct(u, t);
  transform.A[3][1] = dotProduct(v, t);
  transform.A[3][2] = dotProduct(w, t);
  transform.A[3][3] = 1;

  return transform;
}

Matrix
Camera::DeviceTransform(int width, int height)
{
  Matrix transform;

  transform.A[0][0] = transform.A[3][0] = width/2;
  transform.A[1][1] = transform.A[3][1] = height/2;
  transform.A[2][2] = transform.A[3][3] = 1;

  return transform;
  // Transform each vertex of the triangle individually
  // For an n x m image:
  // Image Space coordinate (x, y, z) -> Device Space coordinate (x', y', z'), where:
  // x' = n*(x+1)/2 = nx/2 + n/2
  // y' = m*(y+1)/2 = my/2 + m/2
  // z' = z
  //(This is just scaling [-1, 1] range to [0, screen boundary] range)
  // (x, y, z, 1) * n/2 0   0   0
  //                0   m/2 0   0
  //                0   0   1   0
  //                n/2 m/2 0   1
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

double CalculatePhongShading(LightingParameters& lp, double* viewDirection, double* normal)
{
  double x = dotProduct(lp.lightDir, normal);
  double diffuse = abs(x);

  double R[3] = {0};
  double V[3] = {0};
  scaleVector(2*x, normal, R);
  vectorDifference(R, lp.lightDir, R);
  normalize(R, R);
  normalize(viewDirection, V);
  double d = dotProduct(R, V);
  double specular = 0;
  if (d >= 0)
    specular = pow(d, lp.alpha);

  double phong = lp.Ka + lp.Kd*diffuse + lp.Ks*specular;
  //double phong = lp.Ka;
  //double phong = lp.Kd * diffuse;
  //double phong = lp.Ks * specular;
  return phong;
}

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
    return c;
}


int main()
{
   int width = 1000;
   int height = 1000;

   vtkImageData *image = NewImage(width, height);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = width*height;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.width = width;
   screen.height = height;

   double zBuffer[npixels];

   LightingParameters lp;

   //for (int i = 0; i < 1000; ++i)
   for (int i = 0; i < 1000; i+=1000)
   {
     for (int k = 0; k < npixels*3; ++k)
         buffer[k] = 0;
     for (int k = 0; k < npixels; ++k)
         zBuffer[k] = -1;

     Camera c = GetCamera(i, 1000);
     char filename[64];
     sprintf(filename, "frame%03d", i);

     for (Triangle t : triangles)
     {
       double viewDirection[3] = {0};
       for (int i = 0; i < 3; ++i)
       {
         viewDirection[0] = c.position[0] - t.X[i];
         viewDirection[1] = c.position[1] - t.Y[i];
         viewDirection[2] = c.position[2] - t.Z[i];
         t.shading[i] = CalculatePhongShading(lp, viewDirection, t.normals[i]);
       }

       Matrix viewTransform = c.ViewTransform();
       Matrix cameraTransform = c.CameraTransform();
       Matrix deviceTransform = c.DeviceTransform(width, height);

       // Compose transformations into a single matrix
       Matrix fullTransform = Matrix::ComposeMatrices(cameraTransform, viewTransform);
       fullTransform = Matrix::ComposeMatrices(fullTransform, deviceTransform);

       Triangle tImage = t;
       for (int i = 0; i <= 2; ++i)
       {
         double coordVector[4] = {tImage.X[i], tImage.Y[i], tImage.Z[i], 1};
         double imageVector[4] = {0};

         fullTransform.TransformPoint(coordVector, imageVector);

         double w = imageVector[3];
         tImage.X[i] = imageVector[0]/w;
         tImage.Y[i] = imageVector[1]/w;
         tImage.Z[i] = imageVector[2]/w;
       }

       // Sort verticies by Y-coordinate
       tImage.ySortVertices();

       // Possible way to clean up code: move rasterization code to a separate function, which takes t as a param (so you don't have to refactor t. -> tImage.)

       // Generalized scanline algorithm for triangles
       // The line segments x_0,x_1 and x_0,x_2 are used to calc X-bounds until row == Y[1], after which point x_1,x_2 and x_0,x_2 are used
       // These geometric assumptions can be made because the vertices are sorted by Y-coordinate
       for (int i = 0; i < 2; ++i)
       {
         // Establish max and min rows, accounting for screen boundaries
         int rowMin = ceil_441(tImage.Y[i]);
         if (rowMin < 0)
           rowMin = 0;
         int rowMax = floor_441(tImage.Y[i+1]);
         if (rowMax >= screen.height)
           rowMax = screen.height-1;

         for (int row = rowMin; row <= rowMax; ++row)  // TODO: move unnecessary step reps out of this loop
         {
           double bound1X = 0;
           double bound2X = 0;
           double bound1Z = 0;
           double bound2Z = 0;
           double bound1Color[3] = {0};
           double bound2Color[3] = {0};
           double bound1Shading = 0;
           double bound2Shading = 0;

           // Interpolate x, z, and rgb values for scanline bound 1 on the line segment x_0 x_2,
           // which spans the height of the triangle and will always be one of the sides bounding the scanline
           if (tImage.X[0] == tImage.X[2]) // Account for vertical lines
           {
             bound1X = tImage.X[0];
           }
           // Note: because the vertices are sorted, it is geometrically impossible for Y[0]==Y[2] to be true,
           // so we don't have to check (does not hold for a "triangle" whose vertices are all along one line)
           else
           {
             bound1X = interpolate(tImage.Y[0], tImage.X[0], tImage.Y[2], tImage.X[2], row); // row is implicitly cast to double
           }
           bound1Z = interpolate(tImage.Y[0], tImage.Z[0], tImage.Y[2], tImage.Z[2], row);
           bound1Shading = interpolate(tImage.Y[0], tImage.shading[0], tImage.Y[2], tImage.shading[2], row);
           for(int rgbIndex = 0; rgbIndex <= 2; ++rgbIndex)
             bound1Color[rgbIndex] = interpolate(tImage.Y[0], tImage.colors[0][rgbIndex], tImage.Y[2], tImage.colors[2][rgbIndex], row);

           // Calculate xBound2 using the line segment x_0,x_1 row <= Y[1], and the line segment x_1,x_2 if row > Y[1]
           if (tImage.X[i] == tImage.X[i+1])
           {
             bound2X = tImage.X[i];
             bound2Z = interpolate(tImage.Y[i], tImage.Z[i], tImage.Y[i+1], tImage.Z[i+1], row);
             bound2Shading = interpolate(tImage.Y[i], tImage.shading[i], tImage.Y[i+1], tImage.shading[i+1], row);
             for(int rgbIndex = 0; rgbIndex <= 2; ++rgbIndex)
               bound2Color[rgbIndex] = interpolate(tImage.Y[i], tImage.colors[i][rgbIndex], tImage.Y[i+1], tImage.colors[i+1][rgbIndex], row);
           }
           else if (tImage.Y[i] == tImage.Y[i+1]) // Account for horizontal lines
           {
             // bound1X will always lie along x_0,x_2, so in the case that Y[0]==Y[1] or Y[1]==Y[2], vertex 1 will mark the other x-bound of the horizontal line
             bound2X = tImage.X[1];
             bound2Z = interpolate(tImage.X[i], tImage.Z[i], tImage.X[i+1], tImage.Z[i+1], bound2X); // If dY == 0, use X as the input variable
             bound2Shading = interpolate(tImage.X[i], tImage.shading[i], tImage.X[i+1], tImage.shading[i+1], row);
             for(int rgbIndex = 0; rgbIndex <= 2; ++rgbIndex)
               bound2Color[rgbIndex] = interpolate(tImage.X[i], tImage.colors[i][rgbIndex], tImage.X[i+1], tImage.colors[i+1][rgbIndex], bound2X);
           }
           else
           {
             bound2X = interpolate(tImage.Y[i], tImage.X[i], tImage.Y[i+1], tImage.X[i+1], row);
             bound2Z = interpolate(tImage.Y[i], tImage.Z[i], tImage.Y[i+1], tImage.Z[i+1], row);
             bound2Shading = interpolate(tImage.Y[i], tImage.shading[i], tImage.Y[i+1], tImage.shading[i+1], row);
             for(int rgbIndex = 0; rgbIndex <= 2; ++rgbIndex)
               bound2Color[rgbIndex] = interpolate(tImage.Y[i], tImage.colors[i][rgbIndex], tImage.Y[i+1], tImage.colors[i+1][rgbIndex], row);
           }

           int colMin = 0;
           int colMax = 0;
           if (bound1X <= bound2X)
           {
             colMin = ceil_441(bound1X);
             colMax = floor_441(bound2X);
           }
           else
           {
             colMin = ceil_441(bound2X);
             colMax = floor_441(bound1X);
           }
           if (colMin < 0)
             colMin = 0;
           if (colMax >= screen.width)
             colMax = screen.width-1;

           for (int col = colMin; col <= colMax; ++col)
           {
             double pixelZ = interpolate(bound1X, bound1Z, bound2X, bound2Z, col);
             if (pixelZ > zBuffer[row*screen.width + col])
             {
               zBuffer[row*screen.width + col] = pixelZ;
               for (int rgbIndex = 0; rgbIndex <= 2; ++rgbIndex)
               {
                 double tempColorChannel = interpolate(bound1X, bound1Color[rgbIndex], bound2X, bound2Color[rgbIndex], col); // col is implicitly cast to double
                 double pixelShading = interpolate(bound1X, bound1Shading, bound2X, bound2Shading, col);
                 tempColorChannel *= pixelShading;
                 if (tempColorChannel > 1)
                   tempColorChannel = 1;
                 screen.buffer[3*row*screen.width + 3*col + rgbIndex] = ceil_441(255*tempColorChannel);
               }
             }
           }
         }
       }
     }
     WriteImage(image, filename);
   }
}
