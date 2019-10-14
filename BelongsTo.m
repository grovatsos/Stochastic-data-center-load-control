function Belongs = BelongsTo(A,hyper_plane) %checks if hyper_plane is a row of the matrix A
    
    Belongs = 0;
    matrix_size = size(A);
    no_rows = matrix_size(1);
    for i=1:1:no_rows
        if A(i,:) == hyper_plane
            Belongs = 1;
        end
    end

end