function index = getIndexInGrid( grid, value )
% Finds the index of the closest element to a value in a 2 dimensional grid.
% The grid can be horizontal or vertical
% For example: grid = [1 2 3; 1 2 3; 1 2 3], value = 1.8 => index =2
% Another example: grid = [1 1 1; 2 2 2; 3 3 3], value = 4.1 => index = 3
grid_vec = [];

if grid(1,1) == grid(1,2)   %vertical grid
    grid_vec = grid(:,1);
elseif grid(1,1) == grid(2,1)
    grid_vec = grid(1, :);
else
    error('The grid input is invalid');
        
end
   [~, index] = min(abs(grid_vec-value));

end

