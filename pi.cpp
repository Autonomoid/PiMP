#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include <omp.h>

using namespace std;;

// Function prototypes.
double mc(int samples);
double num_int(int steps);
double pythagoras(double x1, double y1, double x2, double y2);

/////////////////////////////////////////////////////////////////////

// Entry-point.
int main(int argc, char** argv)
{
  // Check console arguments.
  if(argc < 5)
  {
      cout << "USAGE:" << endl;
      cout << "Monte-Carlo: pi <threads> mc <samples> <filename>" << endl;
      cout << "Integration: pi <threads> num_int <steps> <filename>" << endl;
      return 1;
  }

  // Get command-line arguments.
  int threads = atoi(argv[1]);
  string method = string(argv[2]);
  int samples = atoi(argv[3]);
  int steps = atoi(argv[3]);
  char* filename = argv[4];

  // Validate number of threads.
  if(threads < 1)
  {
      cout << "Error: number of threads must be greater than zero." << endl;
      return 1;
  }

  // Validate steps
  if(steps < 1)
  {
      cout << "Error: steps must be greater than zero." << endl;
      return 1;
  }

  // Open a file to save results in.
  ofstream results_file(filename);

  // Note the time taken.
  double t = 0.0;

  // Loop from 1 to n threads.
  for (int i=1; i<threads+1; i++)
  {
    // Set the number of threads.
    omp_set_num_threads(i);

    // Run specified method.
    if(method == "mc")
      t = mc(samples);

    else if(method == "num_int")
      t = num_int(steps);

    // Print error message if method is invalid.
    else
    {
      cout << "Error: possible methods are 'mc'' and 'num_int'" << endl;
      return 1;
    }

    // Print the results
    cout.precision(4);
    cout << "threads = " << i << ", time = " << t << endl;

    // Save results to file
    results_file << i << "," << t << endl;
  }

  // Close file.
  results_file.close();

  return 0;
}

/////////////////////////////////////////////////////////////////////

// Return a random double in the range 'min' to 'max'.
double random(double min, double max)
{
  // Seed the random-number generator using the current time.
  //srand(time(NULL));

  double r = (double)rand() / RAND_MAX;
  return min + r * (max - min);
}

/////////////////////////////////////////////////////////////////////

// Return the Euclidean distance between two (x1, y1) and (x2, y2).
double pythagoras(double x1, double y1, double x2, double y2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  return sqrt(dx*dx + dy*dy);
}

/////////////////////////////////////////////////////////////////////

// Caluclate pi using a Monte-Carlo approximation.
double mc(int samples)
{
  double x, y, d, counted_inside = 0.0;

  // Loop over the number of samples.
  // This is where parallelization is potentially beneficial,
  // as each sample is independent of the rest.

  // Start timing.
  double start = omp_get_wtime();

  #pragma omp parallel for reduction(+:counted_inside)
  for(int i=0; i<samples; i++)
  {
    // Select a random point (x, y) from within the square
    // that encloses a unit circle centred at the origin.
    double x = random(-1.0, 1.0);
    double y = random(-1.0, 1.0);

    // Calculate the distance from the point to the origin.
    double d = pythagoras(x, y, 0.0, 0.0);

    // Is the point within the unit circle?
    if(d < 1.0)
      counted_inside += 1.0;
  }

  // Use the ratio of points inside the unit circle / square
  // and the ratio of the area of the circle / square
  // in order to approximate pi.
  double pi = 4.0 * (counted_inside / samples);

  // Stop timing.
  double end = omp_get_wtime();

  // Print result
  cout.precision(15);
  cout << "pi = " << pi << endl;

  // Return the time taken.
  return (end - start);
}

/////////////////////////////////////////////////////////////////////

// Calculate pi using a numerical integration approximation.
// I.e. integrate (1/1+x^2) using the trapeziod rule.
double num_int(int steps)
{
  double total_area, area, x0, x1, y0, y1 = 0.0;

  // Integrate from 0 to 1.
  double upper_limit = 1.0;
  double lower_limit = 0.0;

  // Calculate width of each interval.
  double stepsize = (upper_limit - lower_limit) / steps;

  // Start timing.
  double start = omp_get_wtime();

  // Loop over all steps, adding the area of each trapezoid to the
  // total area under the curve.
  // This is where paraellization has potential to increase speed
  // as each trapezoid can be worked on independently.

  #pragma omp parallel for private(x0, x1, y0, y1, area)
  //#pragma omp parallel for private(x0, x1, y0, y1, area) reduction(+:total_area)
  for(int i=lower_limit; i<steps; i++)
  {
    // Find the interval's boundaries.
    x0 = i * stepsize;
    x1 = x0 + stepsize;

    // Evaluate y at each boundary.
    y0 = 1/(1 + (x0*x0));
    y1 = 1/(1 + (x1*x1));

    // Calculate are of trapezoid.
    area = 0.5 * (y0 + y1) * stepsize;

    total_area += area;
  }

  // Compute pi.
  double pi = 4.0 * total_area;

  // Stop timing.
  double end = omp_get_wtime();

  // Print result.
  cout.precision(15);
  cout << "pi = " << pi << endl;

  // Return the time taken.
  return (end - start);}
