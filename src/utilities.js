function updateDisc(canvas, disc)
{
    context = canvas.getContext('2d');

    center = disc.getCenter();
    radius = disc.getRadius();

    let path = new Path2D();
    path.arc(center.get(0), center.get(1), radius, 0, 2 * Math.PI, false);

    context.fillStyle = '#e37222'; // TUM Orange
    context.fill(path);

    context.lineWidth = 3;
    context.strokeStyle = '#0065bd'; // TUM Blue
    context.stroke(path);

    return;
}

function callbackFactory(canvas, disc)
{
    context = canvas.getContext('2d');
    points = []

    const f = (event) => {
        const canvasRectangle = canvas.getBoundingClientRect();
        const x = event.clientX - canvasRectangle.left;
        const y = event.clientY - canvasRectangle.top;

        console.log("Include ", [x, y]);
        disc.include(x, y);
        points.push([x, y])
        console.log("Disc center: ", disc.getCenter().get(0), disc.getCenter().get(1));
        console.log("Disc radius: ", disc.getRadius());

        context.clearRect(0, 0, canvas.width, canvas.height);
        updateDisc(canvas, disc);
        points.forEach((point, index) => {
            context.fillStyle = '#000000';
            context.fillRect(point[0]-4, point[1]-4, 8, 8);
        });
    }

    return f;
}
