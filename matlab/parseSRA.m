function [] = parseSRA(input_type)

    global WEBJAR3D RUN_DIR;

    if input_type == 1
        dotbracket  = 'data/sra_input.txt';
        prefix = 'sra_original';
    else
        error('Specify input_type');
    end
    
    WEBJAR3D   = '/Users/anton/Dropbox/BGSURNA/Motifs';
    RUN_DIR    = '/Users/anton/Dropbox/BGSURNA/Motifs/Sequences';
    
    ofn = sprintf('loops_%s.csv', prefix);
    fid = fopen(ofn, 'w');

    F = fastaread(dotbracket);    
        
    % get il and hl positions
    [il,hl] = aExtractSSFromDotBracket(strrep(F(1).Sequence,'-','.'));
    
    S = [];
    for j = 2:length(F)
        S(end+1,1:length(F(j).Sequence)) = F(j).Sequence;
    end
    S = char(S);

    asterisks(1:length(F)-1,1) = '*';    
    c = 1:length(F(1).Sequence);

    % process ils    
    for j = 1:length(il(:,1))

        loop = S(:,[il(j,1):il(j,2) il(j,3):il(j,4)]);
        leftStrand  = S(:,c(il(j,1)):c(il(j,2)));
        rightStrand = S(:,c(il(j,3)):c(il(j,4)));
        loop = [leftStrand asterisks rightStrand];

        il_variants = get_sequence_variants(loop);          

        % location, e.g. 15_20_50_55
        loc = sprintf('%i_%i_%i_%i',il(j,1),il(j,2),il(j,3),il(j,4));
        % id like str1_15_20_50_55
        id  = sprintf('%s%i_%s',prefix,i,loc);

        if ~isempty(il_variants)
            output_csv(il_variants, 'il');
        else
            keyboard;
        end
    end

    % process hls
    for j = 1:length(hl(:,1))
        hairpin = S(:,c(hl(j,1)):c(hl(j,2)));
        hl_variants = get_sequence_variants(hairpin);

        loc = sprintf('%i_%i', hl(j,1),hl(j,2));
        id = sprintf('%s%i_%s',prefix,i,loc);

        if ~isempty(hl_variants)
            output_csv(hl_variants, 'hl');
        else
            keyboard;
        end
    end           
    
    fclose(fid);
    fprintf('Done\n');
    

    function [] = output_csv(variants, loop_type)
        for k = 1:length(variants(:,1))
            fprintf(fid, '"%s","%s","%s","%s","%i"\n', prefix, ... 
                         loop_type, loc, variants{k,1}, variants{k,2});
        end                
    end

end

function [result] = get_sequence_variants(S)

    N = length(S(:,1));
    v = cell(1,N);

    for i = 1:N
        v{i} = S(i,:);        
    end

    [a,b,c] = unique(v);
    
    % x contains counts of each loop from a
    [x,y] = histc(c,1:length(a));
    
    result = {};
    
    for i = 1:length(a)
        seq = a{i};
        seq = strrep(seq,'T','U');
        seq = strrep(seq,'-','');
        
        if ~isempty(strfind(seq,'*'))
            [strands, q] = regexp(seq,'*','split');
            if length(strands{1}) >=2 && length(strands{2}) >= 2 %both strands at least 2 characters
                if length(seq) > 5 % an IL must be at least 2 by 3 + '*' = 6 characters
                    result{end+1,1} = seq;
                    result{end,2} = x(i);            
                end
            end
        elseif length(seq) >=4 % hairpin
            result{end+1,1} = seq;
            result{end,2} = x(i);
        end
    end
    
end
