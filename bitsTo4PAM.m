
function X = bitsTo4PAM (bits)
  %% 00 ->  3
  %% 01 ->  1
  %% 11 -> -1
  %% 10 -> -3
  X = zeros(length(bits)/2  , 1);
  index = 1;
  for i=1 : 2 : length(bits)
    if( bits(i) == 0 && bits(i+1) == 0 )
        X(index) = 3; 
        
    elseif( bits(i) == 0 && bits(i+1) == 1 )  
        X(index) = 1;
         
    elseif( bits(i) == 1 && bits(i+1) == 1 )  
        X(index) = -1;
       
    else  
       X(index) = -3;
    
    endif
  
    index = index + 1;
    
  endfor
  
 
end