function test_suite = test_addFigureTextbox

test_functions = localfunctions();

initTestSuite;

function test_adds_textbox_without_changing_current_axes
    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));
    ax = axes('Parent',fig);
    set(fig,'CurrentAxes',ax);

    textbox = matRad_addFigureTextbox(fig,'RI = 0.9876', ...
        'tag','unitMetricTextbox', ...
        'fontSize',12, ...
        'fontWeight','normal', ...
        'edgeColor',[0.1 0.2 0.3], ...
        'backgroundColor',[0.8 0.9 1.0], ...
        'margin',4);

    assertTrue(ishandle(textbox));
    assertEqual(get(textbox,'Type'),'text');
    assertEqual(get(textbox,'String'),'RI = 0.9876');
    assertEqual(get(textbox,'Tag'),'unitMetricTextbox');
    assertEqual(get(textbox,'FontSize'),12);
    assertEqual(get(textbox,'FontWeight'),'normal');
    assertElementsAlmostEqual(get(textbox,'EdgeColor'),[0.1 0.2 0.3]);
    assertElementsAlmostEqual(get(textbox,'BackgroundColor'),[0.8 0.9 1.0]);
    assertEqual(get(textbox,'Margin'),4);
    assertTrue(get(fig,'CurrentAxes') == ax);

function test_textbox_persists_in_saved_figure
    if exist('savefig','file') ~= 2 || exist('openfig','file') ~= 2
        return;
    end

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));
    matRad_addFigureTextbox(fig,'Gamma-PR = 95.0%', ...
        'tag','savedMetricTextbox');

    figFile = [tempname '.fig'];
    fileCleanup = onCleanup(@() deleteIfExists(figFile));
    savefig(fig,figFile);

    reopenedFig = openfig(figFile,'invisible');
    reopenedCleanup = onCleanup(@() close(reopenedFig));
    textbox = findall(reopenedFig,'Tag','savedMetricTextbox');

    assertEqual(numel(textbox),1);
    assertEqual(get(textbox,'Type'),'text');
    assertEqual(get(textbox,'String'),'Gamma-PR = 95.0%');

function deleteIfExists(filePath)
    if exist(filePath,'file') == 2
        delete(filePath);
    end
