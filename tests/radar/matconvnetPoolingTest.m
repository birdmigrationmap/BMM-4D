classdef matconvnetPoolingTest < matlab.unittest.TestCase
    % MatConvNet's max pooling (used by MistNet) must keep negative maxima.
    % wsrlib maca64 builds with -ffinite-math-only returned 0 instead.

    methods (Test)
        function testMaxPoolingOfNegativeValues(testCase)
            testCase.assumeTrue(exist("vl_nnpool", "file") > 0, "MatConvNet (wsrlib) not on the path.")
            x = -single(reshape(1:16, 4, 4));

            y = vl_nnpool(x, [2 2], 'stride', 2);

            testCase.verifyEqual(y, single([-1 -9; -3 -11]))
        end
    end
end
