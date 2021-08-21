% Threshold ttest results based upon some criteria (e.g., 10 data point in a row.
% 
% inputs:
% p_val
% t_val
% alphalevel -- threshold for deternining initial significance
% numconsectime -- number of consecutive points that must exceed the alphalevel
%   for a point to be retained as significant.
% numconsecchan -- number of consecutive channels that must exceed the alphalevel
% Shlomit Beker 2017 (following a version from the CNL)

function [P T] = SCP_THRESHOLD(p_val,t_val,alphalevel,numconsectime,numconsecchan)


T = zeros(size(t_val));
P = alphalevel*1.3*ones(size(p_val));

for i=1:size(p_val,1)
    
    for k=1:size(p_val,2)
        consec_ch=length(find(p_val(:,1:k) < alphalevel));
        if k>floor(numconsectime/2) && k+(numconsectime/2)<=size(p_val,2)  
            if p_val(i,k-floor(numconsectime/2)+[1:numconsectime]) < alphalevel & ( consec_ch>numconsecchan)
                T(i,k-floor(numconsectime/2)+[1:numconsectime])=t_val(i,k-floor(numconsectime/2)+[1:numconsectime]);
                P(i,k-floor(numconsectime/2)+[1:numconsectime])=p_val(i,k-floor(numconsectime/2)+[1:numconsectime]);
            end
        end
    end
end


%%
% for i=1:size(p_val,1)-numconsecchan
%     for k=1:size(p_val,2)-numconsectime
%         number = find(p_val(i:i+numconsecchan,k:k+numconsectime) < alphalevel);
%         if number == numel(p_val(i:i+numconsecchan,k:k+numconsectime))
%             T(i:i+numconsecchan,k:k+numconsectime)=t_val(i:i+numconsecchan,k:k+numconsectime);
%             P(i:i+numconsecchan,k:k+numconsectime)=p_val(i:i+numconsecchan,k:k+numconsectime);
%         end
%     end
% end

end
