% [sensors]=GenerateData(100,[500 500],[50 50],20,20,'sensors_data.txt')
% Input:
%       N   : Total sensors
%       Area: [Width Lenght]
%       Tile: [TileOfWidth TileOfLenght]
%       K   : Number of underground sensors
%       CZ  : Code Zero - Calculate for Height (1-50)
%       file: File for saving data
% Output: 
%       sensors: Array of sensors' data
function [sensors]=GenerateData_v1(N,Area,Tile,K,CZ,file)
    rowArea=Area(1);
    colArea=Area(2);
    rowTile=Tile(1);
    colTile=Tile(2);
    num_row = rowArea/rowTile
    num_col = colArea/colTile
    
    sensors = [];
    total_sensor_in_cell = round(N/(num_row*num_col),0)
    if total_sensor_in_cell<1
        fprintf('Wrong Tile so total sensor in a cell equa 0\n');
        return;
    end
    each_rand_col_UG = round( (num_row*num_col) / K , 0);
    count_of_UG_in_eachcol = round(num_col / each_rand_col_UG,0);
    %sensors(1,1) = 'x';
    %sensors(1,2) = 'y';
    %sensors(1,3) = 'h';
    count = 1;
    countOfUG = 0;
    for i=1:num_row 
        rand_col = randi([1 num_col],count_of_UG_in_eachcol,1);
        for j=1:num_col
            rand_x = randi([1+rowTile*(i-1) rowTile*i],total_sensor_in_cell,1);
            rand_y = randi([1+colTile*(j-1) colTile*j],total_sensor_in_cell,1);
            if (ismember(j, rand_col) && countOfUG < K)
                rand_h = randi([1 CZ-1],total_sensor_in_cell,1);
                countOfUG = countOfUG + 1;
            else
                rand_h = randi([CZ 50],total_sensor_in_cell,1);
            end
            for k=1:total_sensor_in_cell
                sensors(count,1) = rand_x(k);
                sensors(count,2) = rand_y(k);
                %sensors(count,3) = rand_h(k);
                
                count = count+1;
            end
            
        end
    end
    
    sensors(randperm(size(sensors, 2)), : );
    sensors(:, 3) = 25;
    p = randperm(N,K)
    for i = 1: K
        sensors(p(i), 3) = 20;
    end
    csvwrite(file,sensors);
end

