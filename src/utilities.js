function updateDisc(canvas, disc)
{
    context = canvas.getContext('2d');

    const offset = (point) => [
        point.get(0) + canvas.width / 2,
        canvas.height / 2 - point.get(1)
    ];

    center = offset(disc.getCenter());
    radius = disc.getRadius();

    let path = new Path2D();
    path.arc(center[0], center[1], radius, 0, 2 * Math.PI, false);

    context.fillStyle = '#e37222'; // TUM Orange
    context.fill(path);

    context.lineWidth = 3;
    context.strokeStyle = '#0065bd'; // TUM Blue
    context.stroke(path);

    return;
}

function callbackFactory(canvas, disc)
{
    const top = canvas.offsetTop + canvas.clientTop;
    const left = canvas.offsetLeft + canvas.clientLeft;
    context = canvas.getContext('2d');

    const f = (event) => {
        const x = event.pageX - left / 2;
        const y = top - event.pageY / 2;

        console.log("Include ", [x, y]);
        disc.include(x, y);
        console.log("Disc center: ", disc.getCenter().get(0), disc.getCenter().get(1));
        console.log("Disc radius: ", disc.getRadius());
        updateDisc(canvas, disc);
    }

    return f;
}