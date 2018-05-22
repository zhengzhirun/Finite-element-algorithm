function [g,t] = cvt_2d_sampling ( g_num, it_num, s_num )
%
%  Initialize the generators.
%
  g = rand ( g_num, 2 ) .* 2 - 1;
%
%  Carry out the iteration.
%
  step = 1 : it_num;
  e = nan ( it_num, 1 );
  gm = nan ( it_num, 1 );

  for it = 1 : it_num
%
%  Compute the Delaunay triangle information T for the current nodes.
%
    t = delaunay ( g(:,1), g(:,2) );
%
%  Display the Voronoi cells.
%
    subplot ( 2, 2, 1 );    
    voronoi ( g(:,1), g(:,2), t );
    title_string = sprintf ( 'Voronoi, step %d', it );
    title ( title_string );
    axis equal
    axis square
    drawnow
%
%  Display the Delaunay triangulation.
%
    subplot ( 2, 2, 2 );
    trimesh ( t, g(:,1), g(:,2), zeros(g_num,1), 'EdgeColor', 'r' )
    hold on
    scatter ( g(:,1), g(:,2), 'k.' )
    title_string = sprintf ( 'Delaunay, step %d', it );
    title ( title_string );
    axis square
    view ( 2 )
    hold off
%
%  Generate sample points.
%  New option for fixed grid sampling.
%
    if ( 0 )
      s = rand ( s_num, 2 );
    else
      s2 = floor ( sqrt ( s_num ) );
      s2 = max ( s2, 2 );
      s_num = s2 * s2;
      [ sx, sy ] = meshgrid ( linspace ( 0.0, 1.0, s2 ) );
      sx = reshape ( sx, s2 * s2, 1 );
      sy = reshape ( sy, s2 * s2, 1 );
      s = [ sx, sy ];
    end
%
%  For each sample point, find K, the index of the nearest generator.
%  We do this efficiently by using the Delaunay information with
%  Matlab's DSEARCH command, rather than a brute force nearest neighbor
%  computation.  Also, DSEARCH has been removed, so we need DSEARCHN.
%  
    k = dsearchn ( g, t, s );

    m(:,1) = accumarray ( k, ones(s_num,1) );

    g_new(:,1) = accumarray ( k, s(:,1) ) ./ m(:,1);
    g_new(:,2) = accumarray ( k, s(:,2) ) ./ m(:,1);
%
%  Compute the energy.
%
    e(it,1) = sum ( ( s(:,1) - g(k(:,1),1) ).^2 ...
                  + ( s(:,2) - g(k(:,1),2) ).^2 ) / s_num;
%
%  Display the energy.
%
    subplot ( 2, 2, 3 )
    plot ( step, log ( e ), 'm-*' )
    title ( 'Log (Energy)' )
    xlabel ( 'Step' )
    ylabel ( 'Energy' )
    grid
%
%  Compute the generator motion.
%
    gm(it,1) = sum ( ( g_new(:,1) - g(:,1) ).^2 ...
                   + ( g_new(:,2) - g(:,2) ).^2 ) / g_num;
%
%  Display the generator motion.
%
    subplot ( 2, 2, 4 )
    plot ( step, log ( gm ), 'm-*' )
    title ( 'Log (Average generator motion)' )
    xlabel ( 'Step' )
    ylabel ( 'Energy' )
    grid
%
%  Update the generators.
%
    g = g_new;

  end
end