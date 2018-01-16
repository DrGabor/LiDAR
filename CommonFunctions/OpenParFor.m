function [ flag ] = OpenParFor()
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)   % if exists no parallel computing handle, create one.
    poolsize = parpool;
    flag = 0;
    % disp('Create a Parallel Computing Pool to accelarate process' );
else
    poolsize = p.NumWorkers;
    flag = 1;
    % disp(sprintf( 'Parpool has %d workers', poolsize) );
end
end

