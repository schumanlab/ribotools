classdef BedLine < handle
    % Bed12 format line parser
    
    properties
        chrom
        chromStart
        chromEnd
        name
        transcript
        gene
        score
        strand
        thickStart
        thickEnd
        itemRgb
        blocks
        blockStarts
        blockSizes
        exonStart
        exonEnd
        offset
        cdsStart
        cdsEnd
        cdsSpan
        geneSpan
        linearCoverage
    end
    
    methods
        
        function obj = BedLine()
            
            
        end
        
        function obj = parse(obj,line)
            txt = textscan(line, '%s %n %n %s %n %c %n %n %s %n %s %s','delimiter', '\t');
            obj.chrom = char(txt{1});
            obj.chromStart = txt{2};
            obj.chromEnd = txt{3};
            obj.name = char(txt{4});
            tmp = textscan(obj.name, '%s %s', 'delimiter', ';');
            obj.transcript = char(tmp{1});
            obj.gene = char(tmp{2});
            obj.score = txt{5};
            obj.strand = char(txt{6});
            obj.thickStart = txt{7};
            obj.thickEnd = txt{8};
            obj.itemRgb = char(txt{9});
            obj.blocks = txt{10};
            obj.blockSizes = sscanf(char(txt{11}),'%d,');
            obj.blockStarts = sscanf(char(txt{12}),'%d,');
            obj.exonStart = obj.chromStart + obj.blockStarts;
            obj.exonEnd = obj.exonStart + obj.blockSizes;
            obj.offset = cumsum([0;obj.blockSizes(1:end-1)]);
            obj.geneSpan = sum(obj.blockSizes);
            idxStart = (obj.exonStart <= obj.thickStart) & (obj.thickStart <= obj.exonEnd);
            idxEnd = (obj.exonStart <= obj.thickEnd) & (obj.thickEnd <= obj.exonEnd);
            obj.cdsStart = obj.thickStart - obj.exonStart(idxStart) + obj.offset(idxStart);
            obj.cdsEnd = obj.thickEnd - obj.exonStart(idxEnd) + obj.offset(idxEnd);
            if (obj.strand == '-')
                
                tmpStart = obj.geneSpan - obj.cdsEnd;
                tmpEnd = obj.geneSpan - obj.cdsStart;
                
                obj.cdsStart = tmpStart;
                obj.cdsEnd = tmpEnd;
                
            end
            
            obj.cdsSpan = obj.cdsEnd - obj.cdsStart;
            
            obj.linearCoverage = zeros(obj.geneSpan, 1);
        end
        
        function obj = find(obj, gbed)
    
            for k = 1 : obj.blocks

                idxUse = (gbed.chromStart <= obj.exonEnd(k)) & (obj.exonStart(k) <= gbed.chromEnd);

                if sum(idxUse) == 0
                    continue;
                end
                
                isectStart = max(gbed.chromStart(idxUse), repmat(obj.exonStart(k), sum(idxUse), 1));
                isectEnd = min(gbed.chromEnd(idxUse), repmat(obj.exonEnd(k), sum(idxUse), 1));

                linearStart = isectStart - obj.exonStart(k) + obj.offset(k);
                linearEnd = isectEnd - obj.exonStart(k) + obj.offset(k);
                linearDepth = gbed.depth(idxUse);
                
                
                for j = 1 : size(linearDepth, 1)
                    obj.linearCoverage(linearStart(j) + 1 : linearEnd(j)) = linearDepth(j);
                end
                
            end
            
            % flip for strand
            if (obj.strand == '-')
                obj.linearCoverage = flipud(obj.linearCoverage);
                tmp = obj.cdsStart;
                obj.cdsStart = obj.geneSpan - obj.cdsEnd;
                obj.cdsEnd = obj.geneSpan - tmp;
            end
          
        end
        
        
        function obj = addTranscriptomeTrack(obj, gbed)
            
            for k = 1 : length(gbed.depth)
                obj.linearCoverage(gbed.chromStart(k)+1:gbed.chromEnd(k)) = gbed.depth(k);
            end
            
        end
        
    end
end

