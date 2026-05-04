function matRad_compareRegressionStructs(expected,actual,path,tolerance)
% Recursive comparison with NaNs treated as equal for deterministic outputs.

if nargin < 3 || isempty(path)
    path = 'result';
end
if nargin < 4 || isempty(tolerance)
    tolerance = 1e-12;
end

compareValues(expected,actual,path,tolerance);

end

function compareValues(expected,actual,path,tolerance)
    if isstruct(expected)
        assertTrue(isstruct(actual),sprintf('%s: actual value is not a struct',path));
        assertEqual(size(expected),size(actual),sprintf('%s: struct size differs',path));
        expectedFields = fieldnames(expected);
        actualFields = fieldnames(actual);
        assertEqual(sort(expectedFields),sort(actualFields), ...
            sprintf('%s: struct fields differ',path));
        for elementIx = 1:numel(expected)
            for fieldIx = 1:numel(expectedFields)
                fieldName = expectedFields{fieldIx};
                compareValues(expected(elementIx).(fieldName), ...
                    actual(elementIx).(fieldName), ...
                    sprintf('%s.%s',path,fieldName),tolerance);
            end
        end
    elseif iscell(expected)
        assertTrue(iscell(actual),sprintf('%s: actual value is not a cell',path));
        assertEqual(size(expected),size(actual),sprintf('%s: cell size differs',path));
        for i = 1:numel(expected)
            compareValues(expected{i},actual{i},sprintf('%s{%d}',path,i),tolerance);
        end
    elseif isnumeric(expected) || islogical(expected)
        compareNumeric(expected,actual,path,tolerance);
    elseif ischar(expected)
        assertTrue(ischar(actual),sprintf('%s: actual value is not a char',path));
        assertEqual(expected,actual,sprintf('%s: char values differ',path));
    else
        assertEqual(expected,actual,sprintf('%s: values differ',path));
    end
end

function compareNumeric(expected,actual,path,tolerance)
    assertTrue(isnumeric(actual) || islogical(actual), ...
        sprintf('%s: actual value is not numeric/logical',path));
    assertEqual(size(expected),size(actual),sprintf('%s: numeric size differs',path));

    expected = full(expected);
    actual = full(actual);

    if islogical(expected)
        assertEqual(expected,logical(actual),sprintf('%s: logical values differ',path));
        return;
    end

    expectedNaN = isnan(expected);
    actualNaN = isnan(actual);
    assertEqual(expectedNaN,actualNaN,sprintf('%s: NaN mask differs',path));

    compareIx = ~expectedNaN;
    if any(compareIx(:))
        assertElementsAlmostEqual(actual(compareIx),expected(compareIx), ...
            'absolute',tolerance,sprintf('%s: numeric values differ',path));
    end
end
