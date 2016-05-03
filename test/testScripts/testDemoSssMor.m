classdef testDemoSssMor < sssTest
    
    methods (Test)  
        function mainFunctionality(testCase)
            Opts.pause=false;
            Opts.test=true;
            sss_gettingStarted(Opts); % in sssMor
            sssMOR_gettingStarted(Opts);
        end
    end
end