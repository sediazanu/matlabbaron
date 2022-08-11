classdef loadderivs
    methods (Static)

    function F = loadderivs1(y)
        % define the non-dimensional parameters
        % and give them global values so they can be edited
     global a b c d e
     %input y is a vector holding the prey and predator species 
     dudt = a*y(1) +e*(y(1))^2-b*y(1)*y(2);
     dvdt = c*y(1)*y(2)-d*y(2)*y(2);

     F = [dudt; dvdt];

    end
    
    function F = loadderivs2(t,y)
      % define the non-dimensional parameters
        % and give them global values so they can be edited
    global a b c d e
      %input y is a vector holding the prey and predator species, t
      %represents time
     dudt = a*y(1) +e*(y(1))^2-b*y(1)*y(2);
     dvdt = c*y(1)*y(2)-d*y(2)*y(2);

     F = [dudt; dvdt];

    end
    
    function F = loadderivs3(t,y,Z)
      % define the non-dimensional parameters
        % and give them global values so they can be edited
     global a b c d e
     
      %input y is a vector holding the prey and predator species, t
      %represents time. Z represents the vector holding the increase in
      %delay time.
     delay = Z(:,1);
     dudt = a*y(1) +e*(y(1))^2-b*y(1)*y(2)*delay(1);
     dvdt = c*y(1)*y(2)-d*y(2)*y(2);

     F = [dudt; dvdt];

end
    
    
    end
    
end