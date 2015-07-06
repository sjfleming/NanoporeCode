classdef PoreParams
    %POREPARAMS Static functions for loading and saving parameter files
    %           It would make sense to have the params themselves be
    %           classes, except that it makes passing them to mex-files
    %           slightly more complicated and less universal, so I did it
    %           this way instead.
    
    methods (Static)
        
        function params = Default()
            % Return a params struct with default values
            
            params = [];
            % general aligner parameters
            params.stripe_width = 150;
            params.lik_offset = 4.5;
            % these are the actual probabilities that go into the model
            params.insert_t = 0.02;
            params.insert_c = 0.02;
            params.skip_t = 0.06;
            params.skip_c = 0.06;
            params.stay_t = 0.10;
            params.stay_c = 0.10;
            params.extend_t = 0.14;
            params.extend_c = 0.14;

        end
        
        function params = Load(filename)
            % Loads params from a text file with above fields
            
            params = PoreParams.Default();
            
            if ~exist(filename,'file')
                error(['Parameter file ' filename ' not found']);
            end
            
            fid = fopen(filename);
            tline = fgetl(fid);
            while ischar(tline)
                % parse the file for 'name=val' tokens, with spaces on
                % either side of = sign
                tokens = regexp(tline,'(\w+)\s*=\s*([-+]?[0-9]*\.?[0-9]+)','tokens');
                params.(tokens{1}{1}) = str2double(tokens{1}{2});
                tline = fgetl(fid);
            end
            fclose(fid);
            
        end
        
        function Save(filename, params)
            % Saves params to a text file
            
            fid = fopen(filename,'w');
            fields = fieldnames(params);
            for i=1:numel(fields)
                fprintf(fid,'%s = %f\n',fields{i},params.(fields{i}));
            end
            fclose(fid);
        end
        
    end
end
%}