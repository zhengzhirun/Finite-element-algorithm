function p = cvt_square_nonuniform ( n, sample_num )

%*****************************************************************************80
%
%% CVT_SQUARE_NONUNIFORM computes a CVT in a square, with nonuniform density.
%
%  Discussion:
%
%    For a given density rho(x,y) over the square [-1,+1]x[-1,+1],
%    carry out an iteration for the weighted Centroidal Voronoi 
%    Tessellation (CVT).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    06 November 2016
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Qiang Du, Vance Faber, Max Gunzburger,
%    Centroidal Voronoi Tessellations: Applications and Algorithms,
%    SIAM Review, 
%    Volume 41, 1999, pages 637-676.
%
%  Parameters:
%
%    Input, integer N, the number of generators.
%
%    Input, integer SAMPLE_NUM, the number of sample points.
%
%    Output, real P(N,2), the location of the generators.
%
  timestamp ( )
  fprintf ( 1, '\n' );
  fprintf ( 1, 'CVT_SQUARE_NONUNIFORM:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  A simple demonstration of a CVT computation\n' );
  fprintf ( 1, '  (Centroidal Voronoi Tessellation)\n' );
  fprintf ( 1, '  in a square, with a nonuniform density.\n' );
  
  if ( nargin < 1 )
    n = 100;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CVT_SQUARE_NONUNIFORM - Note:\n' );
    fprintf ( 1, '  No value of N was supplied.\n' );
    fprintf ( 1, '  N is the number of generators.\n' );
    fprintf ( 1, '  A default value N = %d will be used.\n', n );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  User specified number of generators = %d\n',  n );
  end

  if ( nargin < 2 )
    sample_num = 1000 * n;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CVT_SQUARE_NONUNIFORM - Note:\n' );
    fprintf ( 1, '  No value of SAMPLE_NUM was supplied.\n' );
    fprintf ( 1, '  SAMPLE_NUM is the number of sample points.\n' );
    fprintf ( 1, '  A default value SAMPLE_NUM = %d will be used.\n', ...
      sample_num );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  User specified number of sample points = %d\n', ...
      sample_num );
  end
%
%  Clear the figure screen, if already open.
%
  clf
%
%  Randomize the generators in the square [-1,+1]x[-1,+1].
%
  p = square_uniform ( n );

  plot ( p(:,1), p(:,2), 'b.' );
  axis ( [ -1.05, 1.05, -1.05, 1.05 ] )
  line ( [ -1.0, 1.0, 1.0, -1.0, -1.0 ], [ -1.0, -1.0, 1.0, 1.0, -1.0 ], ...
    'Color', 'r', 'LineWidth', 2 );
  title_string = sprintf ( 'Initial Generators' );
  title ( title_string );
  axis equal
  drawnow
  
  for it = 0 : 100
%
%  Compute the Delaunay triangle information T for the current nodes.
%
    t = delaunay ( p(:,1), p(:,2) );
%
%  Display the Voronoi diagram.
%
    voronoi ( p(:,1), p(:,2), t );
    axis ( [ -1.05, 1.05, -1.05, 1.05 ] )
    line ( [ -1.0, 1.0, 1.0, -1.0, -1.0 ], [ -1.0, -1.0, 1.0, 1.0, -1.0 ], ...
      'Color', 'r', 'LineWidth', 2 );
    title_string = sprintf ( 'Voronoi, step %d', it );
    title ( title_string );
    axis equal
    axis tight
    drawnow
%
%  Generate sample points uniformly.  
%    
    s = square_uniform ( sample_num );
%
%  Get the density of each sample point.
%
    d = square_density ( sample_num, s );
%
%  DSEARCHN finds K, the index of the nearest generator for each sample.
% 
    k = dsearchn ( p, t, s );
%
%  Estimate the mass and the center of mass.
%
    mass(1:n,1) = 0.0;
    com(1:n,1) = 0.0;
    com(1:n,2) = 0.0;

    for i = 1 : sample_num
      j = k(i);
      mass(j,1) = mass(j,1) + d(i);
      com(j,1) = com(j,1) + d(i) * s(i,1);
      com(j,2) = com(j,2) + d(i) * s(i,2);
    end

    for i = 1 : n
      if ( mass(i,1) ~= 0.0 )
        com(i,1) = com(i,1) ./ mass(i,1);
        com(i,2) = com(i,2) ./ mass(i,1);
      else
        com(i,1) = p(i,1);
        com(i,2) = p(i,2);
      end
    end
%
%  Replace generators by centers of mass.
%
    p(1:n,1:2) = com(1:n,1:2);
   
  end
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'CVT_SQUARE_NONUNIFORM:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( )

  return
end
function d = square_density ( n, xy )

%*****************************************************************************80
%
%% SQUARE_DENSITY evaluates the density function in the unit square.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    06 June 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of points.
%
%    Input, real XY(N,2), the evaluation points.
%
%    Output, real D(N,1), the density values.
%
  d(1:n,1) = ( xy(1:n,1).^4 + xy(1:n,2).^4 );

  return
end
function p = square_uniform ( n )

%*****************************************************************************80
%
%% SQUARE_UNIFORM returns sample points from the square [-1,+1]x[-1,+1].
%
%  Discussion:
%
%    This routine returns N points sampled uniformly at random
%    from within the unit square.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    06 June 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of points to generate.
%
%    Output, real P(N,2), the sample points.
%
  p = 2.0 * rand ( n, 2 ) - 1.0;

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
