function nodeNamed(graph, name) {
    const node = graph.nodes().filter(node => node.data('name') === name);
    expect(node, `node ${name}`).to.have.length(1);
    return node;
}

function setSimulationParameters() {
    cy.get('#SimParametersForm').should('be.visible');
    cy.get('#numcycles').type('5');
    cy.get('#numtimepts').type('5');
    cy.get('#submitSimParamButton').click();
}

function addNode(name, type, x, y, boundaryCondition) {
    cy.get('#node-type').select(type);
    if (boundaryCondition) {
        cy.get('#boundary-condition-type').select(boundaryCondition);
    }
    cy.get('#node-name').clear().type(name);
    cy.get('#cy').click(x, y);
    cy.window().should(win => {
        expect(nodeNamed(win.cy, name).data('type')).to.equal(type);
    });
}

function clickNode(name) {
    cy.get('#cy').scrollIntoView();
    cy.window().then(win => {
        const container = win.cy.container();
        const position = nodeNamed(win.cy, name).renderedPosition();
        cy.wrap(container).click(
            container.clientLeft + position.x,
            container.clientTop + position.y,
            { scrollBehavior: false }
        );
    });
}

function addVessel(name, x, y) {
    addNode(name, 'vessel', x, y);
    clickNode(name);
    cy.get('#vesselForm').should('be.visible');
    cy.get('#vesselLengthInput').clear().type('1');
    cy.get('#vesselRadiusInput').clear().type('1');
    cy.get('#vesselStenosisDiameterInput').clear().type('0');
    cy.get('#submitVesselButton').click();
}

function mouseAtNode(name, event, buttons) {
    cy.window().then(win => {
        const container = win.cy.container();
        const rect = container.getBoundingClientRect();
        const position = nodeNamed(win.cy, name).renderedPosition();
        cy.wrap(container).trigger(event, {
            eventConstructor: 'MouseEvent',
            clientX: rect.left + container.clientLeft + position.x,
            clientY: rect.top + container.clientTop + position.y,
            button: 0,
            buttons,
            which: 1,
            scrollBehavior: false,
        });
    });
}

function connectNodes(source, target) {
    cy.get('#draw-on').click();
    cy.get('#cy').scrollIntoView();
    mouseAtNode(source, 'mousedown', 1);
    cy.window().should(win => {
        expect(nodeNamed(win.cy, source).hasClass('eh-source'), `start at ${source}`).to.be.true;
    });
    mouseAtNode(target, 'mousemove', 1);
    // Release only after edgehandles has selected this target. Moving a DOM
    // node is not a success condition for a gesture that draws an edge.
    cy.window().should(win => {
        const edges = nodeNamed(win.cy, source).edgesTo(nodeNamed(win.cy, target));
        expect(edges.filter('.eh-preview'), `preview ${source} -> ${target}`).to.have.length(1);
    });
    mouseAtNode(target, 'mouseup', 0);
    cy.window().should(win => {
        const edges = nodeNamed(win.cy, source).edgesTo(nodeNamed(win.cy, target));
        expect(edges, `connection ${source} -> ${target}`).to.have.length(1);
        expect(edges.hasClass('eh-preview'), 'connection is committed').to.be.false;
    });
}

function expectEdges(expected) {
    cy.window().should(win => {
        const actual = win.cy.edges().map(edge => [
            edge.source().data('name'), edge.target().data('name'),
        ]);
        expect(actual, 'directed graph connections').to.have.deep.members(expected);
    });
}

function exportGraph() {
    cy.get('#export-json').click();
    cy.get('@alert').should('not.have.been.called');
    // Observe the actual download boundary; the export handler calls a local
    // function, so stubbing window.downloadJSON does not intercept the export.
    return cy.get('@exportBlob').should('have.been.calledOnce')
        .then(spy => spy.firstCall.args[0].text())
        .then(JSON.parse);
}

beforeEach(() => {
    cy.visit('/');
    cy.window().then(win => {
        cy.stub(win, 'alert').as('alert');
        cy.spy(win.URL, 'createObjectURL').as('exportBlob');
    });
});

describe('Simple Render Tests', () => {
    it('home page rendered', () => {
        cy.get('#cy').should('be.visible');
        cy.get('.collapsible').should('be.visible');
        cy.get('#SimParametersForm').should('be.visible');
    });

    it('Check that Different Nodes are Showing', () => {
        setSimulationParameters();
        addNode('INFLOW', 'boundary_condition', 80, 100, 'FLOW');
        addVessel('vessel0', 200, 100);
        cy.get('.draggable').should('have.length', 2);
        cy.get('@alert').should('not.have.been.called');
    });
});

describe('Node interaction', () => {
    beforeEach(setSimulationParameters);

    it('One Edge Creation', () => {
        addNode('INFLOW', 'boundary_condition', 80, 100, 'FLOW');
        addVessel('vessel0', 200, 100);
        connectNodes('INFLOW', 'vessel0');
        expectEdges([['INFLOW', 'vessel0']]);
        cy.get('@alert').should('not.have.been.called');
    });

    it('Correct Alert was Raised', () => {
        cy.get('#export-json').click();
        cy.get('@alert').should('have.been.calledOnceWithExactly',
            'The model needs at least two boundary conditions');
        cy.get('@exportBlob').should('not.have.been.called');
    });

    for (const [width, height] of [[1000, 660], [1920, 1080]]) {
        it(`Inflow -> Vessel -> OUT (${width}x${height})`, () => {
            cy.viewport(width, height);
            addNode('INFLOW', 'boundary_condition', 80, 100, 'FLOW');
            addVessel('vessel0', 200, 100);
            connectNodes('INFLOW', 'vessel0');
            // Keep the outlet close to the vessel, as in the original failure.
            // Changing node type also changes the controls above the canvas.
            addNode('OUT', 'boundary_condition', 220, 100, 'RESISTANCE');
            connectNodes('vessel0', 'OUT');
            expectEdges([['INFLOW', 'vessel0'], ['vessel0', 'OUT']]);
            exportGraph().then(data => {
                expect(data.vessels).to.have.length(1);
                expect(data.vessels[0].boundary_conditions).to.deep.equal({ inlet: 'INFLOW', outlet: 'OUT' });
                expect(data.boundary_conditions.map(bc => bc.bc_name)).to.deep.equal(['INFLOW', 'OUT']);
            });
        });
    }

    it('Inflow -> Vessel -> Junction -> Vessel -> OUT', () => {
        addNode('INFLOW', 'boundary_condition', 80, 100, 'FLOW');
        addVessel('vessel0', 200, 100);
        connectNodes('INFLOW', 'vessel0');
        addNode('J0', 'junction', 320, 100);
        clickNode('J0');
        cy.get('#junctionForm').should('be.visible');
        cy.get('#submitJunctionButton').click();
        connectNodes('vessel0', 'J0');
        addVessel('vessel1', 440, 100);
        connectNodes('J0', 'vessel1');
        addNode('OUT', 'boundary_condition', 560, 100, 'RESISTANCE');
        connectNodes('vessel1', 'OUT');
        expectEdges([
            ['INFLOW', 'vessel0'], ['vessel0', 'J0'],
            ['J0', 'vessel1'], ['vessel1', 'OUT'],
        ]);
        exportGraph().then(data => {
            expect(data.vessels).to.have.length(2);
            expect(data.vessels[0].boundary_conditions).to.deep.equal({ inlet: 'INFLOW' });
            expect(data.vessels[1].boundary_conditions).to.deep.equal({ outlet: 'OUT' });
            expect(data.junctions).to.deep.equal([{
                junction_name: 'J0', junction_type: 'NORMAL_JUNCTION',
                inlet_vessels: [0], outlet_vessels: [1],
            }]);
        });
    });

    it('Rejects export when the outlet connection is missing', () => {
        addNode('INFLOW', 'boundary_condition', 80, 100, 'FLOW');
        addVessel('vessel0', 200, 100);
        connectNodes('INFLOW', 'vessel0');
        addNode('OUT', 'boundary_condition', 320, 100, 'RESISTANCE');
        expectEdges([['INFLOW', 'vessel0']]);
        cy.get('#export-json').click();
        cy.get('@alert').should('have.been.calledOnceWithExactly',
            'Vessel vessel0 does not have exactly two connections\n' +
            'Boundary condition OUT does not have exactly one connection');
        cy.get('@exportBlob').should('not.have.been.called');
    });
});
