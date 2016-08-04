    function positions = hwayToPos(headways)
        positions = zeros(size(headways));
        for ih = 2:numCars
            positions(ih) = positions(ih - 1) + headways(ih);
        end
        positions = circshift(positions, 1);
    end

