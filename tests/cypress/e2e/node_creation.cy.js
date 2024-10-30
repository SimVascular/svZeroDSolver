import '@4tw/cypress-drag-drop';

describe('Simple Render Tests', () => {
    it('home page rendered', () => {
        cy.visit('/');
        cy.get('#cy', { timeout: 5000 });
        cy.get('.collapsible', { timeout: 5000 });
    });

    it('Check that Different Nodes are Showing', () => {
        cy.visit('/');
        // Enter Simulation Parameters
        cy.get('#SimParametersForm').should('be.visible');
        cy.get('#numcycles').type('5');
        cy.get('#numtimepts').type('5');
        cy.get('#submitSimParamButton').click();

        // Create a BC
        cy.get("#node-type").select("boundary_condition").invoke("val").should("eq", "boundary_condition");
        // Select Inflow
        cy.get("#boundary-condition-type").select("FLOW").invoke("val").should("eq", "FLOW");
        // type node name
        cy.get('#node-name').type('INFLOW');
        // add node
        cy.get('body').click(150, 400);

        // Makes sure the node was added
        cy.get('.draggable').should('have.length', 1);

        cy.get('#node-name').clear();

        // Create a vessel
        cy.get("#node-type").select("vessel").invoke("val").should("eq", "vessel");
        // type node name
        cy.get('#node-name').type('vessel0');
        // add node
        cy.get('body').click(200, 500);
        cy.get('body').click(200, 500);
        cy.get('#vesselForm').should('be.visible');
        cy.get('#vesselLengthInput').type('1');
        cy.get('#vesselRadiusInput').type('1');
        cy.get('#vesselStenosisDiameterInput').type('0');
        cy.get('#submitVesselButton').click();

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 2);

    });
});

describe('Node interaction', () => {
    beforeEach(() => {
        cy.visit('/');
        cy.window().then((win) => {
                cy.stub(win, 'downloadJSON').as('downloadJSONStub');
            });
    });

    it('One Edge Creation', () => {
        let pos1, pos2;

        cy.get('#SimParametersForm').should('be.visible');
        cy.get('#numcycles').type('5');
        cy.get('#numtimepts').type('5');
        cy.get('#submitSimParamButton').click();

        // Create a BC
        cy.get("#node-type").select("boundary_condition").invoke("val").should("eq", "boundary_condition");
        // Select Inflow
        cy.get("#boundary-condition-type").select("FLOW").invoke("val").should("eq", "FLOW");
        // type node name
        cy.get('#node-name').type('INFLOW');

        cy.get('#draw-on').click();

        cy.get('#cy').click(100, 100); // Coordinates for the first node

        cy.get('#node-name').clear();

        // Create a vessel
        cy.get("#node-type").select("vessel").invoke("val").should("eq", "vessel");
        // type node name
        cy.get('#node-name').type('vessel0');
        // add node
        cy.get('body').click(200, 500);
        cy.get('body').click(200, 500);
        cy.get('#vesselForm').should('be.visible');
        cy.get('#vesselLengthInput').type('1');
        cy.get('#vesselRadiusInput').type('1');
        cy.get('#vesselStenosisDiameterInput').type('0');
        cy.get('#submitVesselButton').click();

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 2);



        cy.get('.draggable').first().drag('.draggable:eq(1)', { force: true }).then(() => {
            cy.get('.draggable').eq(1).then($node2 => {
                const rect2 = $node2[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect2.x + rect2.width / 2,
                    clientY: rect2.y + rect2.height / 2,
                    force: true
                });

            });
        });

    });

    it('Correct Alert was Raised', () => {
        // Listen for window alerts
        const alertText = 'The model needs at least two boundary conditions';

        cy.on('window:alert', (str) => {
            // Check if the alert message matches the expected message
            expect(str).to.equal(alertText);
        });

        // Click the export button
        cy.get('#export-json').click();

    });

     it('Inflow -> Vessel -> OUT', () => {

        cy.get('#SimParametersForm').should('be.visible');
        cy.get('#numcycles').type('5');
        cy.get('#numtimepts').type('5');
        cy.get('#submitSimParamButton').click();

        // Create a BC
        cy.get("#node-type").select("boundary_condition").invoke("val").should("eq", "boundary_condition");
        // Select Inflow
        cy.get("#boundary-condition-type").select("FLOW").invoke("val").should("eq", "FLOW");
        // type node name
        cy.get('#node-name').type('INFLOW');

        cy.get('#draw-on').click();

        cy.get('#cy').click(100, 100); // Coordinates for the first node

        cy.get('#node-name').clear();

        // Create a vessel
        cy.get("#node-type").select("vessel").invoke("val").should("eq", "vessel");
        // type node name
        cy.get('#node-name').type('vessel0');
        // add node
        cy.get('body').click(200, 500);
        cy.get('body').click(200, 500);
        cy.get('#vesselForm').should('be.visible');
        cy.get('#vesselLengthInput').type('1');
        cy.get('#vesselRadiusInput').type('1');
        cy.get('#vesselStenosisDiameterInput').type('0');
        cy.get('#submitVesselButton').click();

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 2);

        cy.get('.draggable').first().drag('.draggable:eq(1)', { force: true }).then(() => {
            cy.get('.draggable').eq(1).then($node2 => {
                const rect2 = $node2[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect2.x + rect2.width / 2,
                    clientY: rect2.y + rect2.height / 2,
                    force: true
                });

            });
        });
        cy.get('#node-name').clear();

        // Create a BC
        cy.get("#node-type").select("boundary_condition").invoke("val").should("eq", "boundary_condition");
        // Select Inflow
        cy.get("#boundary-condition-type").select("RESISTANCE").invoke("val").should("eq", "RESISTANCE");
        // type node name
        cy.get('#node-name').type('OUT');

        cy.get('body').click(220, 500);

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 3);

        cy.get('.draggable').eq(1).drag('.draggable:eq(2)', { force: true }).then(() => {
            cy.get('.draggable').eq(2).then($node3 => {
                const rect3 = $node3[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect3.x + rect3.width / 2,
                    clientY: rect3.y + rect3.height / 2,
                    force: true
                });

            });
        });

        // Listen for window alerts and fail the test if any alert is detected
        cy.on('window:alert', (str) => {
            throw new Error(`Unexpected alert: ${str}`);
        });

        // Click the export button
        cy.get('#export-json').click();

    });

    it('Inflow -> Vessel -> Junction -> Vessel -> OUT', () => {

        cy.get('#SimParametersForm').should('be.visible');
        cy.get('#numcycles').type('5');
        cy.get('#numtimepts').type('5');
        cy.get('#submitSimParamButton').click();

        // Create a BC
        cy.get("#node-type").select("boundary_condition").invoke("val").should("eq", "boundary_condition");
        // Select Inflow
        cy.get("#boundary-condition-type").select("FLOW").invoke("val").should("eq", "FLOW");
        // type node name
        cy.get('#node-name').type('INFLOW');

        cy.get('#draw-on').click();

        cy.get('#cy').click(100, 100); // Coordinates for the first node

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 1);

        cy.get('#node-name').clear();

        // Create a vessel
        cy.get("#node-type").select("vessel").invoke("val").should("eq", "vessel");
        // type node name
        cy.get('#node-name').type('vessel0');
        // add node
        cy.get('#cy').click(150, 100);
        cy.get('#cy').click(150, 100);
        cy.get('#vesselForm').should('be.visible');
        cy.get('#vesselLengthInput').type('1');
        cy.get('#vesselRadiusInput').type('1');
        cy.get('#vesselStenosisDiameterInput').type('0');
        cy.get('#submitVesselButton').click();

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 2);

        cy.get('.draggable').first().drag('.draggable:eq(1)', { force: true }).then(() => {
            cy.get('.draggable').eq(1).then($node2 => {
                const rect2 = $node2[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect2.x + rect2.width / 2,
                    clientY: rect2.y + rect2.height / 2,
                    force: true
                });

            });
        });

        cy.get('#node-name').clear();
        // Create a junction
        cy.get("#node-type").select("junction").invoke("val").should("eq", "junction");
        // type node name
        cy.get('#node-name').type('J0');
        cy.get('#cy').click(250, 100)

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 3);

        cy.get('.draggable').eq(1).drag('.draggable:eq(2)', { force: true }).then(() => {
            cy.get('.draggable').eq(2).then($node3 => {
                const rect3 = $node3[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect3.x + rect3.width / 2,
                    clientY: rect3.y + rect3.height / 2,
                    force: true
                });

            });
        });


        cy.get('#node-name').clear();

        // Create a vessel
        cy.get("#node-type").select("vessel").invoke("val").should("eq", "vessel");
        // type node name
        cy.get('#node-name').type('vessel1');
        // add node
        cy.get('#cy').click(300, 100);
        cy.get('#cy').click(300, 100);
        cy.get('#vesselForm').should('be.visible');
        cy.get('#vesselLengthInput').type('1');
        cy.get('#vesselRadiusInput').type('1');
        cy.get('#vesselStenosisDiameterInput').type('0');
        cy.get('#submitVesselButton').click();

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 4);

        cy.get('.draggable').eq(2).drag('.draggable:eq(3)', { force: true }).then(() => {
            cy.get('.draggable').eq(2).then($node3 => {
                const rect3 = $node3[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect3.x + rect3.width / 2,
                    clientY: rect3.y + rect3.height / 2,
                    force: true
                });

            });
        });

        cy.get('#node-name').clear();

        // Create a BC
        cy.get("#node-type").select("boundary_condition").invoke("val").should("eq", "boundary_condition");
        // Select Inflow
        cy.get("#boundary-condition-type").select("RESISTANCE").invoke("val").should("eq", "RESISTANCE");
        // type node name
        cy.get('#node-name').type('OUT');

        cy.get('#cy').click(400, 100);

        // Make sure the node was added
        cy.get('.draggable').should('have.length', 5);

        cy.get('.draggable').eq(3).drag('.draggable:eq(4)', { force: true }).then(() => {
            cy.get('.draggable').eq(2).then($node3 => {
                const rect3 = $node3[0].getBoundingClientRect();
                cy.get('body').trigger('mouseup', {
                    which: 1,
                    clientX: rect3.x + rect3.width / 2,
                    clientY: rect3.y + rect3.height / 2,
                    force: true
                });

            });
        });

        // Listen for window alerts and fail the test if any alert is detected
        cy.on('window:alert', (str) => {
            throw new Error(`Unexpected alert: ${str}`);
        });

        // Click the export button
        cy.get('#export-json').click();

    });
});