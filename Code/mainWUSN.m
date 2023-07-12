% Scenarios:
%   1. S1: 100 nodes, 10% UG, 10% clusters
%   2. S2: 100 nodes, 20% UG, 10% clusters
%   3. S3: 100 nodes, 30% UG, 10% clusters
%   4. S4: 200 nodes, 10% UG, 10% clusters
%   5. S5: 200 nodes, 20% UG, 10% clusters
%   6. S6: 200 nodes, 30% UG, 10% clusters
%
%   X,C,m,Eps,maxTest, BS, alpha, kU

function mainWUSN()
    scenarios = {
                    [100 0.1 0.1 500 500 50 50 20],
                    [100 0.1 0.2 500 500 50 50 20],
                    [100 0.1 0.3 500 500 50 50 20],
                    [200 0.1 0.1 500 500 50 50 20],
                    [200 0.1 0.2 500 500 50 50 20],
                    [200 0.1 0.3 500 500 50 50 20]
                };
    number_of_scenarios = size(scenarios,1);
    number_of_round = 10;
    m = 2;
    Eps = 0.01;
    maxTest = 1;
    BS = [250 250 25];
    alpha = 0.002;
    running_times = [];
    for k = 1:number_of_scenarios(1)
        N = scenarios{k,1}(1);
        C = scenarios{k,1}(2)*N;
        kU = scenarios{k,1}(3)*N;
        Area = [scenarios{k,1}(4) scenarios{k,1}(5)];
        Tile = [scenarios{k,1}(6) scenarios{k,1}(7)];
        CZ = scenarios{k,1}(8);
        inputfile = sprintf('%s\\data\\sensors_data_scenario%d.txt',pwd,k);
        
        GenerateData(N,Area,Tile,kU,CZ,inputfile);
        
        
        X=dlmread(inputfile);
        outputFCMWUSN_EC_rounds = [];
        outputFCM_EC_rounds = [];
        
        outputFCMWUSN_LN_rounds = [];
        outputFCM_LN_rounds = [];
        
        outputFCMWUSN_LUN_rounds = [];
        outputFCM_LUN_rounds = [];
        
        for r = 1: number_of_round    
            starting_time = now;           
            output = FCMChuanWUSN(X,C,m,Eps,maxTest, BS, alpha, kU);
           
            outputFCMWUSN_EC_rounds{r,1} = output{1,1};
            outputFCM_EC_rounds{r,1} = output{2,1};
            outputFCMWUSN_LN_rounds{r,1} = output{3,1};
            outputFCM_LN_rounds{r,1} = output{4,1};
            outputFCMWUSN_LUN_rounds{r,1} = output{5,1};
            outputFCM_LUN_rounds{r,1} = output{6,1};
            
            total_seconds = now - starting_time;
            running_times(k,r) = total_seconds;
        end
        
        outputfile = {
                        sprintf('%s\\output\\sensors_result_FCMWUSN_EC_scenario%d_%s.csv',pwd,k,datestr(now, 'yyyymmdd_HHMMSS')),
                        sprintf('%s\\output\\sensors_result_FCM_EC_scenario%d_%s.csv',pwd,k,datestr(now, 'yyyymmdd_HHMMSS')),
                        sprintf('%s\\output\\sensors_result_FCMWUSN_LN_scenario%d_%s.csv',pwd,k,datestr(now, 'yyyymmdd_HHMMSS')),
                        sprintf('%s\\output\\sensors_result_FCM_LN_scenario%d_%s.csv',pwd,k,datestr(now, 'yyyymmdd_HHMMSS')),
                        sprintf('%s\\output\\sensors_result_FCMWUSN_LUN_scenario%d_%s.csv',pwd,k,datestr(now, 'yyyymmdd_HHMMSS')),
                        sprintf('%s\\output\\sensors_result_FCM_LUN_scenario%d_%s.csv',pwd,k,datestr(now, 'yyyymmdd_HHMMSS'))
                        };
        writecell(outputFCMWUSN_EC_rounds, outputfile{1,1});
        writecell(outputFCM_EC_rounds, outputfile{2,1});
        writecell(outputFCMWUSN_LN_rounds, outputfile{3,1});
        writecell(outputFCM_LN_rounds, outputfile{4,1});
        writecell(outputFCMWUSN_LUN_rounds, outputfile{5,1});
        writecell(outputFCM_LUN_rounds, outputfile{6,1});
    end
    dlmwrite(sprintf('%s\\output\\running_time.txt',pwd),running_times);
end

